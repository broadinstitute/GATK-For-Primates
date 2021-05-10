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

import "./structs/structs.wdl"
import "./tasks/CollectInfo.wdl" as CollectInfo
import "./tasks/FASTQtoBAM.wdl" as FASTQtoBAM
import "./tasks/BaseRecalibrator.wdl" as BQSR
import "./tasks/GatherVcfs.wdl" as Gather
import "./tasks/HardFilter.wdl" as HardFilter
import "./tasks/ProduceFinalCallset.wdl" as ProduceFinalCallset
import "./tasks/VariantCallAndGenotype.wdl" as VariantCallAndGenotype
import "./tasks/VariantRecalibrator.wdl" as VQSR
import "./tasks/GenerateTSVs.wdl" as GenerateTSVs

## WORKFLOW DEFINITION
workflow GATKForPrimatesGermlineSNPsIndels_GATK4 {

    String pipeline_version = "pre-alpha"
    
    ####################################################################
    ## Initial workflow setup and configuration
    ####################################################################

    ## Workflow inputs
    input {

        ## Collect workflow options from input JSON
        String mode # Options: initial repeat final
        Boolean bwamem2 = false # Options: true / false to use either bwa (as bwa mem) or bwamem2 (as bwamem2 mem)
        Boolean flowcell_patterned = false # Options: true / false to use either bwa (as bwa mem) or bwamem2 (as bwamem2 mem)
        Int? merge_contigs_into_num_partitions # Option to speed up GenomicsDBImport with large number of contigs
        File? truth_set_SNPs # Option to provide a SNP truth set for VQSR (training set is produced via hard filtering)
        File? truth_set_INDELs # Option to provide an INDEL truth set for VQSR (training set is produced via hard filtering)

        ## Collect reference input file
        File ref

        ## Define arrays from input JSON; structural definitions are at the end of workflow
        Array[sampleInfo] sampleList
        Array[scatterInfo] scatterList
        
        ## Placeholders for when I add runtime/docker info
        ## Is it just me or has the 'Genomes in the Cloud' docker from Broad not been updated
        ## in three years? https://hub.docker.com/r/broadinstitute/genomes-in-the-cloud
        ## So using this one instead for bwa and samtools
        String docker_image_gatk = "broadinstitute/gatk:4.2.0.0"
        String docker_image_bwa_and_samtools = "bioslimcontainers/bwa-samtools:bwa-0.7.17_samtools-1.11"
        String docker_image_python = "python:latest"

    }

    ## Determine correct reference index files depending on version of bwa is being used
    File ref_dict = sub(ref, ".fa", ".dict")
    Array[File] ref_idxs = if bwamem2 then prefix(ref + ".", ["fai", "amb", "ann", "pac", "0123","bwt.2bit.64"]) else prefix(ref + ".", ["fai", "amb", "ann", "pac", "bwt", "sa"])

    ## Set correct command for executing bwa depending on which version is being used
    String execute_bwa = if bwamem2 then "bwa-mem2" else "bwa"

    ## Set configuration Booleans as true or false
    Boolean input_is_fastq = mode == "initial"
    Boolean truth_SNPs = defined(truth_set_SNPs)
    Boolean truth_INDELs = defined(truth_set_INDELs)

    ####################################################################
    ## Initial scatters to collect variables needed later scattering
    ####################################################################

    ## Scatter over scatters to get list of scatter names and intervals
    ## All this is necessary to sub-scatter later on
    scatter (sc in scatterList) {
        call collectScatters {
            input:
                scatterName = sc.name,
                scatterIntervals = sc.intervals,
        }
    }

    ## Define arrays comprising all scatter names and intervals
    ## Array[String] allScatterNames = collectScatters.allScatterNames
    ## Array[String] allScatterIntervals = collectScatters.allScatterIntervals

    ## Scatter over scatters to get list of group and sample names
    ## All this is necessary to sub-scatter later on and for counting number of samples
    scatter (sample in sampleList) {
        call collectAllGroupNames {
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

        ## Scatter over samples
        scatter (sample in sampleList) {

            ## Map paired-end reads into BAM file
            call FastqToBwaMem {
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
                docker_image = docker_image_bwa_and_samtools,
            }

            ## Run MarkDuplicatesSpark
            ## Outputs BAM file sorted by co-ordinate
            call MarkDuplicatesSpark {
            input:
                sampleName = sample.name,
                input_bam = FastqToBwaMem.output_bam,
                flowcell_patterned = flowcell_patterned,
                docker_image = docker_image_gatk,
            }

            ## Run SetNmMdAndUqTags
            ## Outputs BAM file with fixed tags
            call SortAndFixTags {
            input:
                ref = ref,
                ref_dict = ref_dict,
                ref_idxs = ref_idxs,
                sampleName = sample.name,
                sampleGroup = sample.group,                
                input_bam = MarkDuplicatesSpark.output_bam,
                docker_image = docker_image_gatk,
            }

        }

        ## Generate a TSV (not really a 'map' per se) of the newly created BAM files
        call GenerateMapOfNewBAMs {
        input:
            new_bams = SortAndFixTags.output_bam,
            new_bam_groups = SortAndFixTags.output_group,
            new_bam_names = SortAndFixTags.output_name,
            new_bam_indexes = SortAndFixTags.output_bam_index,
        }

        ## Get array of all group names for scattering over later
        #Array[String] allGroupsForScattering = GenerateMapOfNewBAMs.allGroupsForScattering

        ## Coerce the file from File? into File because Cromwell cannot accommodate
        ## using files that may not exist in ternary operators (i.e. if-then-else)
        call coerceFile as coerceMapOfNewBAMs {
        input:
            input_file = GenerateMapOfNewBAMs.map,
        }

        ## Read the new TSV into a map
        Array[Array[String]] map_of_new_bams = read_tsv(coerceMapOfNewBAMs.coerced)

    }

    ####################################################################
    ## MODE = REPEAT/FINAL: Use existing BAM files
    ####################################################################

    if (!input_is_fastq) {

        ## Collect info on each sample from the input JSON file
        scatter (sample in sampleList) {
            call collectExistingBamInfo {
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
        call GenerateMapOfExistingBAMs {
        input:
            existing_bam_names = existingSamples,
            existing_bam_groups = existingGroups,
            existing_bams = existingBAMs,
            existing_bam_indexes = existingBAIs,
        }

        ## Coerce the file from File? into File because Cromwell cannot accommodate
        ## using files that may not exist in ternary operators (i.e. if-then-else)
        call coerceFile as coerceMapOfExistingBAMs {
        input:
            input_file = GenerateMapOfExistingBAMs.map,
        }

        ## Read the new TSV into a map
        Array[Array[String]] map_of_existing_bams = read_tsv(coerceMapOfExistingBAMs.coerced)

    }

    ####################################################################
    ## MODE = ANY: Work out which BAM files to use hereon out
    ####################################################################

    ## Work out which map of BAM files to use; either new or existing
    ## This has to be passed to a task because Cromwell cannot accommodate
    ## using files that may not exist in ternary operators (i.e. if-then-else)
    call FinalizeInputs {
    input:
        map_of_new_bams = map_of_new_bams,
        map_of_existing_bams = map_of_existing_bams,
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
            call HaplotypeCaller {
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
            }

        }
        ## We've just stopped scattering over samples but are still scattering
        ## over the user-defined scatter interval.

        ## Now we scatter over each group as pair.left and send an array of all member samples as pair.right
        scatter (pair in as_pairs(collect_by_key(zip(HaplotypeCaller.output_gvcf_group,HaplotypeCaller.output_gvcf)))) {

            ## Now we import all gVCFS per group into a GenomicsDB
            ## There is a bit of a bottleneck here; because we used pair and collect_by_key we cannot send the
            ## .tbi file with the .gvcf.gz -- so this is created again in the task.
            call GenomicsDBImport {
            input:
                database_name = "initial_variant_calls_" + pair.left + "_" + sc.name,
                input_gvcfs = pair.right,
                scatterIntervals = sc.intervals,
                docker_image = docker_image_gatk,
                merge_contigs_into_num_partitions = merge_contigs_into_num_partitions,
            }

            call GenotypeGenomicsDB {
            input:
                ref = ref,
                ref_dict = ref_dict,
                ref_idxs = ref_idxs,
                input_genomicsdb = GenomicsDBImport.output_genomicsdb,
                database_name = "initial_variant_calls_" + pair.left + "_" + sc.name,
                scatterIntervals = sc.intervals,
                docker_image = docker_image_gatk,
                groupName = pair.left,
            }

            ## This only really needs to be done if it's the final run
            ## or when there is no defined SNP truth set, but WDL
            ## makes this difficult because it doesn't like conditionals
            ## -- and it doesn't take very long anyway
            call HardFilterSNPs {
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
            }

            call HardFilterINDELs {
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
            }

        } ## Stop scattering over group/sample

    } ## Stop scattering over scatter list

    ####################################################################
    ## MODE = INITIAL/REPEAT: Base recalibration
    ####################################################################

    if (mode == "initial" || mode == "repeat") {

        ## Gather hard filtered SNPs
        scatter (pair in as_pairs(collect_by_key(zip(flatten(HardFilterSNPs.output_groupName),flatten(HardFilterSNPs.output_filtered_SNPs))))) {

            call GatherHardFilteredVcfs as GatherSNPs {
            input:
                input_filtered = pair.right,
                ref_dict = ref_dict,
                docker_image = docker_image_gatk,
                groupName = pair.left,
                type = "SNP",
            }
        }

        scatter (pair in as_pairs(collect_by_key(zip(flatten(HardFilterINDELs.output_groupName),flatten(HardFilterINDELs.output_filtered_INDELs))))) {

            call GatherHardFilteredVcfs as GatherINDELs {
            input:
                groupName = pair.left,
                input_filtered = pair.right,
                ref_dict = ref_dict,
                docker_image = docker_image_gatk,
                type = "INDEL",
            }
        }

        ## Then we scatter over each input BAM file to recalibrate
        scatter (col in FinalizeInputs.finalized_inputs) {

            ## Call BaseRecalibrator etc. for each BAM file
            call BaseRecalibrator {
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
            }
        }

    }

    ## Workflow is now complete if mode is initial or repeat
    
    ####################################################################
    ## MODE = FINAL: Variant recalibration/filtration
    ####################################################################

    if (mode == "final" && truth_SNPs) {

        scatter (pair in as_pairs(collect_by_key(zip(flatten(GenotypeGenomicsDB.output_groupName),flatten(GenotypeGenomicsDB.output_genotypes))))) {

            call GatherUnfilteredVcfs as GatherUnfilteredSNPs {
            input:
                groupName = pair.left,
                input_unfiltered = pair.right,
                type = "SNP",
                ref = ref,
                ref_dict = ref_dict,
                ref_idxs = ref_idxs,
                docker_image = docker_image_gatk,
            }

        }

        call VariantRecalibrator as RecalibrateSNPs {
        input:
            ref = ref,
            ref_dict = ref_dict,
            ref_idxs = ref_idxs,
            truth_set = truth_set_SNPs,
            input_sites = GatherUnfilteredSNPs.output_unfiltered_sites_only,
            input_sites_indexes = GatherUnfilteredSNPs.output_unfiltered_sites_only_index,
            type = "SNP",
            docker_image = docker_image_gatk,
        }

    } 

    if (mode == "final" && truth_INDELs) {

        scatter (pair in as_pairs(collect_by_key(zip(flatten(GenotypeGenomicsDB.output_groupName),flatten(GenotypeGenomicsDB.output_genotypes))))) {

            call GatherUnfilteredVcfs as GatherUnfilteredINDELs {
            input:
                groupName = pair.left,
                input_unfiltered = pair.right,
                type = "INDEL",
                ref = ref,
                ref_dict = ref_dict,
                ref_idxs = ref_idxs,
                docker_image = docker_image_gatk,

            }
        } 

        call VariantRecalibrator as RecalibrateINDELs {
        input:
            ref = ref,
            ref_dict = ref_dict,
            ref_idxs = ref_idxs,
            truth_set = truth_set_INDELs,
            input_sites = GatherUnfilteredINDELs.output_unfiltered_sites_only,
            input_sites_indexes = GatherUnfilteredINDELs.output_unfiltered_sites_only_index,
            type = "INDEL",
            docker_image = docker_image_gatk,
        }

    }



####################################################################
   ## MODE = FINAL: Perform final round of variant calls
####################################################################

    call ProduceFinalCallset {
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
    }


    scatter (col in FinalizeInputs.finalized_inputs) {

        ## Call HaplotypeCaller
        ## This calls haplotypes within 100 bp -ip of all the final callset sites
        call HaplotypeCaller as HaplotypeCallerFinal {
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
        }

    }

    call GenomicsDBImportFinal {
    input:
        input_gvcfs = HaplotypeCallerFinal.output_gvcf,
        input_gvcfs_indexes = HaplotypeCallerFinal.output_gvcf_index,
        database_name = "final_variant_calls",
        FinalCallset = ProduceFinalCallset.FinalCallset,
        FinalCallset_index = ProduceFinalCallset.FinalCallset_index,
        docker_image = docker_image_gatk,
    }

    call GenotypeGenomicsDB as GenotypeGenomicsDBFinal {
    input:
        ref = ref,
        ref_dict = ref_dict,
        ref_idxs = ref_idxs,
        input_genomicsdb = GenomicsDBImportFinal.output_genomicsdb,
        database_name = "final_variant_calls",
        FinalCallset = ProduceFinalCallset.FinalCallset,
        FinalCallset_index = ProduceFinalCallset.FinalCallset_index,        
        docker_image = docker_image_gatk,
    }

}
