version development

## Copyright Broad Institute and Wisconsin National Primate Research Center,
## University of Wisconsin-Madison, 2021
## 
## Complete germline short variant discovery pipeline optimized for
## non-human primates, following proposed GATK Best Practices for
## Non-Human Animal Genomes.
##
## Requirements/expectations :
## - Paired-end reads in FASTQ format (for 'initial' mode) or mapped and
##   de-duplicated BAM files with read group information (for 'repeat'
##   and 'final' modes).
## - A reference genome indexed using either bwa or bwa-mem2.
## - An input JSON file listing the reference, the required metadata for
##   each sample, and any optional parameters.
## - Truth SNP and/or INDEL sets, if available.
## - A list of 'scatters' (as described in the paper, DOI TBD).
##
## Outputs :
## - TBD
##
## Software version requirements :
## - bwa 0.7.17
## - GATK 4.2.0.0
## - Samtools 1.12
## - Python 2.7
##
## Cromwell version support 
## - Testing in progress on v58
##
## Runtime parameters are optimized for Broad's Google Cloud Platform implementation.
## For program versions, see docker containers.
##
## LICENSING : 
## This script is released under the WDL source code license (BSD-3) (see LICENSE in 
## https://github.com/broadinstitute/wdl). Note however that the programs it calls may 
## be subject to different licenses. Users are responsible for checking that they are
## authorized to run all programs before running this script. Please see the dockers
## for detailed licensing information pertaining to the included programs.

import "https://raw.githubusercontent.com/broadinstitute/GATK-For-Primates/main/structs/structs.wdl"
import "https://raw.githubusercontent.com/broadinstitute/GATK-For-Primates/main/tasks/CollectInfo.wdl" as CollectInfo
import "https://raw.githubusercontent.com/broadinstitute/GATK-For-Primates/main/tasks/FASTQtoBAM.wdl" as FASTQtoBAM
import "https://raw.githubusercontent.com/broadinstitute/GATK-For-Primates/main/tasks/BaseRecalibrator.wdl" as BQSR
import "https://raw.githubusercontent.com/broadinstitute/GATK-For-Primates/main./tasks/GatherVcfs.wdl" as Gather
import "https://raw.githubusercontent.com/broadinstitute/GATK-For-Primates/main/tasks/HardFilter.wdl" as HardFilter
import "https://raw.githubusercontent.com/broadinstitute/GATK-For-Primates/main/tasks/ProduceFinalCallset.wdl" as ProduceFinalCallset
import "https://raw.githubusercontent.com/broadinstitute/GATK-For-Primates/main/tasks/VariantCallAndGenotype.wdl" as VariantCallAndGenotype
import "https://raw.githubusercontent.com/broadinstitute/GATK-For-Primates/main/tasks/VariantRecalibrator.wdl" as VQSR
import "https://raw.githubusercontent.com/broadinstitute/GATK-For-Primates/main/tasks/GenerateTSVs.wdl" as GenerateTSVs
import "https://raw.githubusercontent.com/broadinstitute/GATK-For-Primates/main/tasks/FinalizeInputs.wdl" as FinalizeInputs

## WORKFLOW DEFINITION
workflow GATKForPrimatesGermlineSNPsIndels_GATK4 {

    String pipeline_version = "pre-alpha"

    input {
        String mode # JSON options: initial / repeat / final
        Boolean bwamem2 = false # JSON options: true / false; indicating bwa (as bwa mem) or bwamem2 (as bwamem2 mem)
        Boolean flowcell_patterned = false # JSON options: true / false
        Int? merge_contigs_into_num_partitions # JSON options: optional parameter for GenomicsDBImport
        File? truth_set_SNPs # JSON options: optional SNP truth set for VQSR (training set is produced via hard filtering)
        File? truth_set_INDELs # JSON options: optional INDEL truth set for VQSR (training set is produced via hard filtering)

        File ref

        ## Define arrays from input JSON; structural definitions are in the structs/structs.wdl file
        Array[sampleInfo] sampleList
        Array[scatterInfo] scatterList
        
        ## Docker containers
        String docker_image_gatk = "broadinstitute/gatk:4.2.0.0"
        String docker_image_bwa_and_samtools = "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.7-1603303710"
        String docker_image_python = "python:latest"
        String docker_image_gatk_with_R_and_ggplot = "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.7-1603303710"

        ## Optional runtime arguments
        Int? preemptible_tries
        File? gatk_override
        String? gatk_docker_override   
    }

    ## Set either bwa or bwamem2 index files and execution command; note bwa-mem2 requires change in docker file
    File ref_dict = sub(ref, ".fa", ".dict")
    Array[File] ref_idxs = if bwamem2 then prefix(ref + ".", ["fai", "amb", "ann", "pac", "0123","bwt.2bit.64"]) else prefix(ref + ".", ["fai", "amb", "ann", "pac", "bwt", "sa"])
    String execute_bwa = if bwamem2 then "bwa-mem2" else "/usr/gitc/bwa"

    ## Set configuration Booleans as true or false
    Boolean input_is_fastq = mode == "initial"
    Boolean truth_SNPs = defined(truth_set_SNPs)
    Boolean truth_INDELs = defined(truth_set_INDELs)

    ####################################################################
    ## Collect scatter, group and sample information from JSON
    ####################################################################

    ## Collect scatter names and intervals
    scatter (sc in scatterList) {
        call CollectInfo.collectScatters as collectScatters {
            input:
                scatterName = sc.name,
                scatterIntervals = sc.intervals,
        }
    }

    ## Collect names of groups and samples within them
    scatter (sample in sampleList) {
        call CollectInfo.collectAllGroupNames as collectAllGroupNames {
        input:
            groupName = sample.group,
            sampleName = sample.name,
        }
    }

    ## Array[String] allGroupsForScattering = collectAllGroupNames.allGroupNames
    Array[String] allSampleNames = collectAllGroupNames.allSampleNames


    ####################################################################
    ## MODE = INITIAL: Map FASTQ files to BAM
    ####################################################################
    
    if (input_is_fastq && mode == "initial") {

        ## Scatter over each sample
        scatter (sample in sampleList) {

            ## Map paired-end reads into BAM file
            call FASTQtoBAM.FastqToBwaMem as FastqToBwaMem {
            input:
                ref = ref,
                ref_dict = ref_dict,
                ref_idxs = ref_idxs,
                sampleName = sample.name,
                rgID = sample.rgID,
                rgSM = sample.rgSM,
                rgLB = sample.rgLB,
                R1 = sample.R1,
                R2 = sample.R2,
                execute_bwa = execute_bwa,
                preemptible_tries = preemptible_tries,
                docker_image = docker_image_bwa_and_samtools,
            }

            ## Run MarkDuplicatesSpark, outputs co-ordinate-sorted BAM
            call FASTQtoBAM.MarkDuplicatesSpark as MarkDuplicatesSpark {
            input:
                sampleName = sample.name,
                input_bam = FastqToBwaMem.output_bam,
                flowcell_patterned = flowcell_patterned,
                preemptible_tries = preemptible_tries,
                docker_image = docker_image_gatk,
            }

            ## Run SetNmMdAndUqTags, outputs BAM with fixed tags
            call FASTQtoBAM.SortAndFixTags as SortAndFixTags {
            input:
                ref = ref,
                ref_dict = ref_dict,
                ref_idxs = ref_idxs,
                sampleName = sample.name,
                sampleGroup = sample.group,                
                input_bam = MarkDuplicatesSpark.output_bam,
                preemptible_tries = preemptible_tries,
                docker_image = docker_image_gatk,
            }

        }

        ## Generate TSV detailing the newly created BAM files
        call GenerateTSVs.BAMs as GenerateTSVOfNewBAMs {
        input:
            bams = SortAndFixTags.output_bam,
            bam_groups = SortAndFixTags.output_group,
            bam_names = SortAndFixTags.output_name,
            bam_indexes = SortAndFixTags.output_bam_index,
            preemptible_tries = preemptible_tries,
        }

        Array[Array[String]] tsv_of_new_bams = read_tsv(GenerateTSVOfNewBAMs.tsv)

    }

    ####################################################################
    ## MODE = REPEAT/FINAL: Use existing BAM files
    ####################################################################

    if (!input_is_fastq) {

        ## Collect info on each sample from the input JSON file
        scatter (sample in sampleList) {
            call CollectInfo.collectExistingBamInfo as collectExistingBamInfo {
                input:
                    sampleName = sample.name,
                    sampleGroup = sample.group,
                    sampleBAM = sample.bam,
                    sampleBAI = sample.bam_index,
            }
        }

        ## Pull all this information together into arrays
        Array[String] existingSamples = collectExistingBamInfo.allNames
        Array[String] existingGroups = collectExistingBamInfo.allGroups
        Array[String] existingBAMs = collectExistingBamInfo.allBams
        Array[String] existingBAIs = collectExistingBamInfo.allBais

        ## Use the arrays to generate a TSV (not really a 'map' per se) of the existing BAM files
        call GenerateTSVs.BAMs as GenerateTSVOfExistingBAMs {
        input:
            bam_names = existingSamples,
            bam_groups = existingGroups,
            bams = existingBAMs,
            bam_indexes = existingBAIs,
        }

        ## Read the new TSV into a map
        Array[Array[String]] tsv_of_existing_bams = read_tsv(GenerateTSVOfExistingBAMs.tsv)

    }

    ####################################################################
    ## MODE = ANY: Work out which BAM files to use hereon out
    ####################################################################

    ## Work out which map of BAM files to use; either new or existing
    ## This has to be passed to a task because Cromwell cannot accommodate
    ## using files that may not exist in ternary operators (i.e. if-then-else)
    call FinalizeInputs.FinalizeInputs as FinalizeInputs {
    input:
        tsv_of_new_bams = tsv_of_new_bams,
        tsv_of_existing_bams = tsv_of_existing_bams,
    }
    
    ####################################################################
    ## MODE = ANY: Run HaplotypeCaller per sample per scatter interval
    ####################################################################

    ## Now we scatter over the user-defined scatter intervals
    scatter (sc in scatterList) {

        ## Then we scatter over each input BAM file
        scatter (col in FinalizeInputs.finalized_inputs) {

            ## Call HaplotypeCaller
            ## This calls haplotypes for the current sample at the current scatter
            call VariantCallAndGenotype.HaplotypeCaller as HaplotypeCaller {
            input:
                ref = ref,
                ref_dict = ref_dict,
                ref_idxs = ref_idxs,
                input_bam = col[2],
                input_bam_index = col[3],
                gvcf_basename = col[0] + ".haplotypes." + col[1] + "." + sc.name,
                sampleName = col[0],
                sampleGroup = col[1],
                scatterName = sc.name,
                scatterIntervals = sub(sc.intervals, " ", " -L "),
                docker_image = docker_image_gatk,
                preemptible_tries = preemptible_tries,
            }

        }
        ## We've just stopped scattering over samples but are still scattering
        ## over the user-defined scatter interval.

        ## Now we scatter over each group as pair.left and send an array of all member samples as pair.right
        scatter (pair in as_pairs(collect_by_key(zip(HaplotypeCaller.output_gvcf_group,HaplotypeCaller.output_gvcf)))) {

            ## Now we import all gVCFS per group into a GenomicsDB
            ## There is a bit of a bottleneck here; because we used pair and collect_by_key we cannot send the
            ## .tbi file with the .gvcf.gz -- so this is created again in the task.
            call VariantCallAndGenotype.GenomicsDBImport as GenomicsDBImport {
            input:
                database_name = "initial_variant_calls_" + pair.left + "_" + sc.name,
                input_gvcfs = pair.right,
                scatterIntervals = sc.intervals,
                docker_image = docker_image_gatk,
                merge_contigs_into_num_partitions = merge_contigs_into_num_partitions,
            }

            call VariantCallAndGenotype.GenotypeGenomicsDB as GenotypeGenomicsDB{
            input:
                ref = ref,
                ref_dict = ref_dict,
                ref_idxs = ref_idxs,
                input_genomicsdb = GenomicsDBImport.output_genomicsdb,
                database_name = "initial_variant_calls_" + pair.left + "_" + sc.name,
                scatterIntervals = sc.intervals,
                groupName = pair.left,
                docker_image = docker_image_gatk,
                preemptible_tries = preemptible_tries,
            }

            ## This only really needs to be done if it's the final run
            ## or when there is no defined SNP truth set, but WDL
            ## makes this difficult because it doesn't like conditionals
            ## -- and it doesn't take very long anyway
            call HardFilter.SNPs as HardFilterSNPs {
            input:
                ref = ref,
                ref_dict = ref_dict,
                ref_idxs = ref_idxs,
                input_genotypes = GenotypeGenomicsDB.output_genotypes,
                input_genotypes_index = GenotypeGenomicsDB.output_genotypes_index,
                groupName = pair.left,
                scatterName = sc.name,
                scatterIntervals = sc.intervals,
                docker_image = docker_image_gatk,
                preemptible_tries = preemptible_tries,
            }

            call HardFilter.INDELs as HardFilterINDELs {
            input:
                ref = ref,
                ref_dict = ref_dict,
                ref_idxs = ref_idxs,
                input_genotypes = GenotypeGenomicsDB.output_genotypes,
                input_genotypes_index = GenotypeGenomicsDB.output_genotypes_index,
                groupName = pair.left,
                scatterName = sc.name,
                scatterIntervals = sc.intervals,
                docker_image = docker_image_gatk,
                preemptible_tries = preemptible_tries,
            }

        } ## Stop scattering over group/sample

    } ## Stop scattering over scatter list

    ####################################################################
    ## MODE = INITIAL/REPEAT: Base recalibration
    ####################################################################

    if (mode == "initial" || mode == "repeat") {

        ## Gather hard filtered SNPs
        scatter (pair in as_pairs(collect_by_key(zip(flatten(HardFilterSNPs.output_groupName),flatten(HardFilterSNPs.output_filtered_SNPs))))) {

            call Gather.HardFilteredVcfs as GatherSNPs {
            input:
                input_filtered = pair.right,
                ref_dict = ref_dict,
                docker_image = docker_image_gatk,
                groupName = pair.left,
                type = "SNP",
                docker_image = docker_image_gatk,
                preemptible_tries = preemptible_tries,
            }
        }

        scatter (pair in as_pairs(collect_by_key(zip(flatten(HardFilterINDELs.output_groupName),flatten(HardFilterINDELs.output_filtered_INDELs))))) {

            call Gather.HardFilteredVcfs as GatherINDELs {
            input:
                groupName = pair.left,
                input_filtered = pair.right,
                ref_dict = ref_dict,
                docker_image = docker_image_gatk,
                type = "INDEL",
                docker_image = docker_image_gatk,
                preemptible_tries = preemptible_tries,
            }
        }

        ## Then we scatter over each input BAM file to recalibrate
        scatter (col in FinalizeInputs.finalized_inputs) {

            ## Call BaseRecalibrator etc. for each BAM file
            call BQSR.BaseRecalibrator as BaseRecalibrator {
            input:
                    ref = ref,
                    ref_dict = ref_dict,
                    ref_idxs = ref_idxs,
                    input_bam = col[2],
                    input_bam_index = col[3],
                    input_SNP_sites = GatherSNPs.output_filtered_sites_only,
                    input_SNP_sites_indexes = GatherSNPs.output_filtered_sites_only_index,
                    input_INDEL_sites = GatherINDELs.output_filtered_sites_only,
                    input_INDEL_sites_indexes = GatherINDELs.output_filtered_sites_only_index,
                    sampleName = col[0],
                    docker_image = docker_image_gatk,
                    preemptible_tries = preemptible_tries,
            }
            
            ## Make plots for each BAM file -- requires separate task as R and ggplot needed
            ## and the genomesinthecloud docker uses an old version of GATK
            call BQSR.AnalyzeCovariates as AnalyzeCovariates {
            input:
                    sampleName = col[0],
                    table_before = BaseRecalibrator.table_before,
                    table_after = BaseRecalibrator.table_after,
                    docker_image = docker_image_gatk_with_R_and_ggplot,
                    preemptible_tries = preemptible_tries,
            }
        }

    }

    ## Workflow is now complete if mode is initial or repeat
    
    ####################################################################
    ## MODE = FINAL: Variant recalibration/filtration
    ####################################################################

    if (mode == "final" && truth_SNPs) {

        scatter (pair in as_pairs(collect_by_key(zip(flatten(GenotypeGenomicsDB.output_groupName),flatten(GenotypeGenomicsDB.output_genotypes))))) {

            call Gather.UnfilteredVcfs as GatherUnfilteredSNPs {
            input:
                groupName = pair.left,
                input_unfiltered = pair.right,
                type = "SNP",
                ref = ref,
                ref_dict = ref_dict,
                ref_idxs = ref_idxs,
                docker_image = docker_image_gatk,
                preemptible_tries = preemptible_tries,
            }

        }

        call VQSR.VariantRecalibrator as RecalibrateSNPs {
        input:
            ref = ref,
            ref_dict = ref_dict,
            ref_idxs = ref_idxs,
            truth_set = truth_set_SNPs,
            input_sites = GatherUnfilteredSNPs.output_unfiltered_sites_only,
            input_sites_indexes = GatherUnfilteredSNPs.output_unfiltered_sites_only_index,
            type = "SNP",
            docker_image = docker_image_gatk,
            preemptible_tries = preemptible_tries,
        }

    } 

    if (mode == "final" && truth_INDELs) {

        scatter (pair in as_pairs(collect_by_key(zip(flatten(GenotypeGenomicsDB.output_groupName),flatten(GenotypeGenomicsDB.output_genotypes))))) {

            call Gather.UnfilteredVcfs as GatherUnfilteredINDELs {
            input:
                groupName = pair.left,
                input_unfiltered = pair.right,
                type = "INDEL",
                ref = ref,
                ref_dict = ref_dict,
                ref_idxs = ref_idxs,
                docker_image = docker_image_gatk,
                preemptible_tries = preemptible_tries,
            }
        } 

        call VQSR.VariantRecalibrator as RecalibrateINDELs {
        input:
            ref = ref,
            ref_dict = ref_dict,
            ref_idxs = ref_idxs,
            truth_set = truth_set_INDELs,
            input_sites = GatherUnfilteredINDELs.output_unfiltered_sites_only,
            input_sites_indexes = GatherUnfilteredINDELs.output_unfiltered_sites_only_index,
            type = "INDEL",
            docker_image = docker_image_gatk,
            preemptible_tries = preemptible_tries,
        }

    }



####################################################################
   ## MODE = FINAL: Perform final round of variant calls
####################################################################

    call ProduceFinalCallset.ProduceFinalCallset as ProduceFinalCallset {
    input:
        ref_dict = ref_dict,
        SNPs_hard_filtered = GatherSNPs.output_filtered_sites_only,
        SNPs_hard_filtered_index = GatherSNPs.output_filtered_sites_only_index,
        INDELs_hard_filtered = GatherINDELs.output_filtered_sites_only,
        INDELs_hard_filtered_index = GatherINDELs.output_filtered_sites_only_index,
        SNPs_recalibrated = RecalibrateSNPs.output_recalibrated_sites_only,
        SNPs_recalibrated_index = RecalibrateSNPs.output_recalibrated_sites_only_index,
        INDELs_recalibrated = RecalibrateINDELs.output_recalibrated_sites_only,
        INDELs_recalibrated_index = RecalibrateINDELs.output_recalibrated_sites_only_index,
        docker_image = docker_image_gatk,
        preemptible_tries = preemptible_tries,
    }

    scatter (col in FinalizeInputs.finalized_inputs) {

        ## Call HaplotypeCaller
        ## This calls haplotypes within 100 bp -ip of all the final callset sites
        call VariantCallAndGenotype.HaplotypeCaller as HaplotypeCallerFinal {
        input:
            ref = ref,
            ref_dict = ref_dict,
            ref_idxs = ref_idxs,
            input_bam = col[2],
            input_bam_index = col[3],
            gvcf_basename = col[0] + ".haplotypes",
            sampleName = col[0],
            FinalCallset = ProduceFinalCallset.FinalCallset,
            FinalCallset_index = ProduceFinalCallset.FinalCallset_index,
            docker_image = docker_image_gatk,
            preemptible_tries = preemptible_tries,
        }

    }

    call VariantCallAndGenotype.GenomicsDBImportFinal as GenomicsDBImportFinal {
    input:
        input_gvcfs = HaplotypeCallerFinal.output_gvcf,
        input_gvcfs_indexes = HaplotypeCallerFinal.output_gvcf_index,
        database_name = "final_variant_calls",
        FinalCallset = ProduceFinalCallset.FinalCallset,
        FinalCallset_index = ProduceFinalCallset.FinalCallset_index,
        docker_image = docker_image_gatk,
        preemptible_tries = preemptible_tries,
    }

    call VariantCallAndGenotype.GenotypeGenomicsDB as GenotypeGenomicsDBFinal {
    input:
        ref = ref,
        ref_dict = ref_dict,
        ref_idxs = ref_idxs,
        input_genomicsdb = GenomicsDBImportFinal.output_genomicsdb,
        database_name = "final_variant_calls",
        FinalCallset = ProduceFinalCallset.FinalCallset,
        FinalCallset_index = ProduceFinalCallset.FinalCallset_index,        
        docker_image = docker_image_gatk,
        preemptible_tries = preemptible_tries,
    }

}
