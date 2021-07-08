version 1.0

## Copyright Broad Institute and Wisconsin National Primate Research Center,
## University of Wisconsin-Madison, 2021
## 
## Complete germline short variant discovery pipeline optimized for
## non-human primates, following proposed GATK Best Practices for
## Non-Human Animal Genomes.
##
## For requirements, expectations and outputs, please review the complete
## documentation at: https://github.com/broadinstitute/GATK-For-Primates
##
## Software version requirements :
## - Cromwell 63
## - bwa 0.7.15 (note from GITC)
## - Samtools 1.11 (note from GITC)
## - GATK 4.2.0.0 (note that GATK 4.1.8.0 is used in GITC)
## - Python 3.9.5
##
## Program versions can be changed by defining alternative containers.
## Runtime parameters are optimized for Terra (https://www.terra.bio/)
##
## LICENSING : 
## This script is released under the WDL source code license (BSD-3) (see LICENSE in 
## https://github.com/broadinstitute/wdl). However, the programs it calls may be
## subject to different licenses. Users are responsible for checking that they are
## authorized to run all programs before running this script. Please see the docker
## containers for detailed licensing information pertaining to the included programs.

import "./structs/structs.wdl"
import "./workflows/ValidateUserInputs.wdl" as preflight
import "./workflows/UnmappedInputToAlignedBAM.wdl" as mapToBam
import "./tasks/Utilities.wdl" as utilities
import "./tasks/HaplotypesToGenotypes.wdl" as bamToVcf
import "./workflows/PipeGenotypes.wdl" as pipeGenotypes
import "./tasks/ProcessVCFs.wdl" as processVCFs
import "./workflows/VariantRecalibration.wdl" as variantRecalibration
import "./workflows/VariantRecalibrationPartial.wdl" as partialVariantRecalibration
import "./workflows/BaseRecalibration.wdl" as baseRecalibration
import "./workflows/functions/CollectByKey.wdl" as collectByKey

##########################################################################
## WORKFLOW DEFINITION
##########################################################################

workflow GATKForPrimatesGermlineSNPsIndels_GATK4 {

    String pipeline_version = "pre-alpha"

    ##########################################################################
    ## DEFINE WORFKLOW INPUTS
    ##########################################################################

    input {

        ## Collect workflow mode from JSON
        String mode # options: initial / repeat / final
 
        ## Collect optional variables from JSON
        Boolean validate_reference_vcf = false # options: true / false; if true will perform ValidateVariants on reference in 'Initial' mode only
        Boolean validate_truth_sets = false # options: true / false; if true will force ValidateVariants on truth sets in 'Repeat' and 'Final' modes
        Boolean bwamem2 = false # options: true / false; indicating bwa (as bwa mem) or bwamem2 (as bwamem2 mem) ***-Coming-Soon-***
        Boolean flowcell_patterned = false # options: true / false; this influences pixel distance when marking duplicates
        Int? merge_contigs_into_num_partitions # options: optional parameter for GenomicsDBImport ***-Coming-Soon-***

        ## Collect reference files from JSON
        File ref
        File ref_dict
        File ref_fai
        File? ref_amb
        File? ref_ann
        File? ref_pac
        File? ref_bwt
        File? ref_sa
        File? ref_0123
        File? ref_bwt_2bit_64

        ## Collect truth sets from JSON
        File? truth_set_SNPs # options: optional SNP truth set (sites-only VCF file) for VQSR (training set is produced via hard filtering)
        File? truth_set_SNPs_index # index for the above
        File? truth_set_INDELs # options: optional INDEL truth set (sites-only VCF file) for VQSR (training set is produced via hard filtering)
        File? truth_set_INDELs_index # index for the above
        
        ## Define arrays from input JSON; definitions are in the structs/structs.wdl file
        Array[sampleInfo]+ sampleList
        Array[scatterInfo]+ scatterList
        
        ## Define containers and paths
        String container_gatk = "broadinstitute/gatk:4.2.0.0"
        String container_gitc = "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.5.7-2021-06-09_16-47-48Z"
        String container_python = "python:3.9.5"

        ## Define paths
        String path_to_gitc = "/usr/gitc/"
        String path_to_gitc_gatk = "/usr/gitc/gatk4/"

    }

    ##########################################################################
    ## Process input options and set default values
    ##########################################################################

    ## Assign correct index files to ref_idxs depending on version of bwa
    Array[File?]+ ref_idxs = if bwamem2 then [ref_amb, ref_ann, ref_pac, ref_0123, ref_bwt_2bit_64] else [ref_amb, ref_ann, ref_pac, ref_bwt, ref_sa]
    ## Set mode Booleans as true or false
    Boolean mode_is_initial = mode == "initial"
    Boolean mode_is_final = mode == "final"
    ## Set Booleans determining which genotype files should be built and piped downstream
    Boolean pipe_unfiltered = if (defined(truth_set_INDELs)) || (defined(truth_set_SNPs)) && (mode_is_final) then true else false
    Boolean pipe_unfiltered_sites_only = if (pipe_unfiltered) then true else false
    Boolean pipe_hard_filtered_SNPs = if ((defined(truth_set_INDELs)) && (mode_is_final)) || ((!defined(truth_set_INDELs)) && (!defined(truth_set_SNPs)) && (mode_is_final)) then true else false
    Boolean pipe_hard_filtered_INDELs = if ((defined(truth_set_SNPs)) && (mode_is_final)) || ((!defined(truth_set_INDELs)) && (!defined(truth_set_SNPs)) && (mode_is_final)) then true else false
    Boolean pipe_hard_filtered_SNPs_sites_only = if ((defined(truth_set_INDELs)) && (!defined(truth_set_SNPs)) && (mode_is_final)) then false else true
    Boolean pipe_hard_filtered_INDELs_sites_only = if ((defined(truth_set_SNPs)) && (!defined(truth_set_INDELs)) && (mode_is_final)) then false else true

    ##########################################################################
    ## Validate user input files and variables
    ##########################################################################

    call preflight.validateUserInputs as validateUserInputs {
        input:
            ref = ref,
            ref_dict = ref_dict,
            ref_fai = ref_fai,
            scatterList = scatterList,
            sampleList = sampleList,
            truth_set_SNPs = truth_set_SNPs,
            truth_set_SNPs_index = truth_set_SNPs_index,
            truth_set_INDELs = truth_set_INDELs,
            truth_set_INDELs_index = truth_set_INDELs_index,
            mode_is_initial = mode_is_initial,
            validate_reference_vcf = validate_reference_vcf,
            validate_truth_sets = validate_truth_sets,
            ## Runtime
            container_gatk = container_gatk,
            container_python = container_python,
    }

    ####################################################################
    ## IF MODE == INITIAL: Map either paired FASTQ or uBAM to BAM
    ####################################################################
    
    if (mode_is_initial) {

        call mapToBam.unmappedInputToAlignedBAM as mapReads {
            input:
                ref = ref,
                ref_dict = ref_dict,
                ref_fai = ref_fai,
                ref_idxs = ref_idxs,
                sampleList = sampleList,
                flowcell_patterned = flowcell_patterned,
                bwamem2 = bwamem2,
                ## Runtime
                container_gatk = container_gatk,
                container_gitc = container_gitc,
                path_to_gitc = path_to_gitc,
                path_to_gitc_gatk = path_to_gitc_gatk,
                container_python = container_python,
        }

    }

    ####################################################################
    ## Generate JSON of BAM files to use downstream
    ####################################################################

    call utilities.generateBamJSON as bamJSON {
        input:
            sampleGroups = select_first([mapReads.sampleGroups,validateUserInputs.groupNames]),
            sampleNames = select_first([mapReads.sampleNames,validateUserInputs.sampleNames]),
            bams = select_first([mapReads.bams,validateUserInputs.bams]),
            bam_indexes = select_first([mapReads.bam_indexes,validateUserInputs.bam_indexes]),
            ## Runtime
            container = container_python,
    }

    Array[bamInfo]+ bamList = read_json(bamJSON.file)

    ####################################################################
    ## Perform germline joint-variant discovery
    ####################################################################

    ## Scatter over each user-defined scatter interval
    scatter (scttr in scatterList) {

        ## Scatter over each BAM file
        scatter (bamFile in bamList) {

            ## Call haplotypes (i.e. per BAM, per scatter)
            call bamToVcf.haplotypeCaller as haplotypeCaller {
                input:
                    ref = ref,
                    ref_dict = ref_dict,
                    ref_fai = ref_fai,
                    input_bam = bamFile.bam,
                    input_bam_index = bamFile.bam_index,
                    sampleGroup = bamFile.sampleGroup,
                    title_gvcf = bamFile.sampleName + "_haplotypes_" + scttr.name,
                    scatterIntervals = sub(scttr.intervals, " ", " -L "),
                    # Runtime
                    container = container_gatk,
            }

        }

        ## Generate haplotype JSON, comprising 'Cohort' if multiple taxa
        call utilities.generateHaplotypeJSON as haplotypeJSON {
            input:
                sampleGroups = haplotypeCaller.output_sampleGroup,
                title_gvcfs = haplotypeCaller.output_title_gvcf,
                gvcfs = haplotypeCaller.output_gvcf,
                gvcf_indexes = haplotypeCaller.output_gvcf_index,
                multiple_taxonomic_groups = validateUserInputs.multiple_taxonomic_groups, 
                # Runtime
                container = container_python,
        }

        ## Read haplotypeJSON as a haplotypeInfo struct and coerce into arrays
        Array[haplotypeInfo]+ haplotypeList = read_json(haplotypeJSON.file)
        scatter (haplos in haplotypeList) {
            String haplotype_groupNames = haplos.groupName
            String haplotype_gvcfs = haplos.gvcf
        }

        ## Scatter over each group, collecting haplotype gVCF files for each
        ## WDL 1.1: scatter (pair in as_pairs(collect_by_key(zip(select_all(haplotype_groupNames),select_all(haplotype_gvcfs))))) {
        ## Code below for WDL 1.0:
        call collectByKey.collectByKey as collectHaplotypesByGroup {
            input:
                allGroups = validateUserInputs.groupNames,
                groups = select_all(haplotype_groupNames),
                members = select_all(haplotype_gvcfs),
        }

        scatter (pair in collectHaplotypesByGroup.collected) {

            ## Restore in WDL 1.1
            #scatter (each_haplotype_gvcf in pair.right) {
            #    String each_haplotype_gvcf_index = each_haplotype_gvcf + ".tbi"
            #}

            scatter (each_haplotype_gvcf in pair.right) {
                String each_haplotype_gvcf_index = each_haplotype_gvcf + ".tbi"
            }
            
            ## Import haplotypes to GenomicsDBs
            call bamToVcf.genomicsDBImport as genomicsDBImport {
                input:
                    input_gvcfs = pair.right,
                    input_gvcfs_indexes = each_haplotype_gvcf_index,
                    title_gdb = "gdb_" + pair.left + "_" + scttr.name,
                    scatterIntervals = scttr.intervals,
                    merge_contigs_into_num_partitions = merge_contigs_into_num_partitions,
                    # Runtime
                    container = container_gatk,
            }

            ## Genotype each database
            call bamToVcf.genotypeGenomicsDB as genotypeGenomicsDB {
                input:
                    ref = ref,
                    ref_dict = ref_dict,
                    ref_fai = ref_fai,
                    input_genomicsdb = genomicsDBImport.output_genomicsdb,
                    input_genomicsdb_title = genomicsDBImport.output_genomicsdb_title,
                    output_vcf_title = "variant_calls_" + pair.left + "_" + scttr.name + "_genotyped",
                    groupName = pair.left,
                    scatterName = scttr.name,
                    scatterIntervals = scttr.intervals,
                    # Runtime
                    container = container_gatk,
            }

        }

    }

    ####################################################################
    ## Produce and pipe the variants that are needed downstream
    ####################################################################

    ## Scatter over each group, collecting genotype VCF files for each
    ## WDL 1.1: scatter (pair in as_pairs(collect_by_key(zip(select_all(flatten(genotypeGenomicsDB.output_groupName)),select_all(flatten(genotypeGenomicsDB.output_genotypes_vcf)))))) {
    ## Code below for WDL 1.0:
    call collectByKey.collectByKey as collectGenotypesByGroup {
        input:
            allGroups = validateUserInputs.groupNames,
            groups = select_all(flatten(genotypeGenomicsDB.output_groupName)),
            members = select_all(flatten(genotypeGenomicsDB.output_genotypes_vcf)),
    }
    Array[Pair[String,Array[String]]] collectedGenotypesByGroup = collectGenotypesByGroup.collected
    
    scatter (pair in collectedGenotypesByGroup) {

        call pipeGenotypes.pipeGenotypes as pipe {
            input:
                ref = ref,
                ref_dict = ref_dict,
                ref_fai = ref_fai,
                gdb_groupName = pair.left,
                gdb_vcfs = pair.right,
                pipe_unfiltered = pipe_unfiltered,
                pipe_unfiltered_sites_only = pipe_unfiltered_sites_only,
                pipe_hard_filtered_SNPs = pipe_hard_filtered_SNPs,
                pipe_hard_filtered_INDELs = pipe_hard_filtered_INDELs,
                pipe_hard_filtered_SNPs_sites_only = pipe_hard_filtered_SNPs_sites_only,
                pipe_hard_filtered_INDELs_sites_only = pipe_hard_filtered_INDELs_sites_only,
                # Runtime
                container = container_gatk,
        }

    }

    call utilities.generateGenotypeJSON as genotypeJSON {
        input:
            groupName = pipe.groupName,
            loc_unfiltered_vcf = pipe.loc_unfiltered_vcf,
            loc_unfiltered_vcf_index  = pipe.loc_unfiltered_vcf_index,
            loc_unfiltered_sites_only_vcf = pipe.loc_unfiltered_sites_only_vcf,
            loc_unfiltered_sites_only_vcf_index = pipe.loc_unfiltered_sites_only_vcf_index,
            loc_hard_filtered_SNPs_vcf = pipe.loc_hard_filtered_SNPs_vcf,
            loc_hard_filtered_SNPs_vcf_index = pipe.loc_hard_filtered_SNPs_vcf_index,
            loc_hard_filtered_INDELs_vcf = pipe.loc_hard_filtered_INDELs_vcf,
            loc_hard_filtered_INDELs_vcf_index = pipe.loc_hard_filtered_INDELs_vcf_index,
            loc_hard_filtered_SNPs_sites_only_vcf = pipe.loc_hard_filtered_SNPs_sites_only_vcf,
            loc_hard_filtered_SNPs_sites_only_vcf_index = pipe.loc_hard_filtered_SNPs_sites_only_vcf_index,
            loc_hard_filtered_INDELs_sites_only_vcf = pipe.loc_hard_filtered_INDELs_sites_only_vcf,
            loc_hard_filtered_INDELs_sites_only_vcf_index = pipe.loc_hard_filtered_INDELs_sites_only_vcf_index,
            # Runtime
            container = container_python,
    }

    Array[genotypeInfo]+ genotypeList = read_json(genotypeJSON.file)

    ####################################################################
    ## IF MODE == INITIAL/REPEAT: PERFORM BASE RECALIBRATION
    ####################################################################

    if (mode_is_initial || !mode_is_final) {    

        ## Scatter over each input BAM file to recalibrate
        scatter (bamFile in bamList) {

            call baseRecalibration.BQSR as BQSR {
                input:
                    ref = ref,
                    ref_dict = ref_dict,
                    ref_fai = ref_fai,
                    input_bam = bamFile.bam,
                    input_bam_index = bamFile.bam_index,
                    input_SNP_sites = select_all(pipe.hard_filtered_SNPs_sites_only_vcf),
                    input_SNP_sites_indexes = select_all(pipe.hard_filtered_SNPs_sites_only_vcf_index),
                    input_INDEL_sites = select_all(pipe.hard_filtered_INDELs_sites_only_vcf),
                    input_INDEL_sites_indexes = select_all(pipe.hard_filtered_INDELs_sites_only_vcf_index),
                    sampleName = bamFile.sampleName,
                    container_gatk = container_gatk,
                    container_gitc = container_gitc,
                    path_to_gitc_gatk = path_to_gitc_gatk,
            }

        }

    }

    ## Workflow is now complete if mode is initial or repeat


    ####################################################################
    ## IF MODE = FINAL: Cannot VQSR; make list of master loci by merging
    ## all sites-only hard-filtered SNPs and INDELs from each group
    ####################################################################

    if ((!defined(truth_set_INDELs)) && (!defined(truth_set_SNPs)) && (mode_is_final)) {

        call processVCFs.mergeManyVCFs as mergeMasterLocifromHF {
            input:
                input_array_1_vcf = select_all(pipe.hard_filtered_SNPs_sites_only_vcf),
                input_array_1_vcf_index = select_all(pipe.hard_filtered_SNPs_sites_only_vcf_index),
                input_array_2_vcf = select_all(pipe.hard_filtered_INDELs_sites_only_vcf),
                input_array_2_vcf_index = select_all(pipe.hard_filtered_INDELs_sites_only_vcf_index),
                ref_dict = ref_dict,
                output_vcf_title = "mergeINDELsFromVQSRwithHardFilteredSNPs",
                # Runtime
                container = container_gatk,
        }
        
    }


    ####################################################################
    ## IF MODE = FINAL: Can VQSR; proceed either fully or partially
    ####################################################################

    ## If mode is final and at least one truth set is available
    if ((defined(truth_set_INDELs)) || (defined(truth_set_SNPs)) && (mode_is_final)) {

        ## Scatter over the available genotype sets per group
        scatter (genotypeFile in genotypeList) {

            ## If both truth sets are available, recalibrate each in series
            if (defined(truth_set_SNPs) && defined(truth_set_INDELs)) {
                call variantRecalibration.VQSR as VQSR {
                    input:
                        ref = ref,
                        ref_dict = ref_dict,
                        ref_fai = ref_fai,
                        groupName = genotypeFile.groupName,
                        unfiltered_vcf = genotypeFile.unfiltered_vcf,
                        unfiltered_vcf_index = genotypeFile.unfiltered_vcf_index,
                        unfiltered_vcf_sites_only = genotypeFile.unfiltered_sites_only_vcf,
                        unfiltered_vcf_sites_only_index = genotypeFile.unfiltered_sites_only_vcf_index,
                        truth_set_INDELs = truth_set_INDELs,
                        truth_set_INDELs_index = truth_set_INDELs_index,
                        truth_set_SNPs = truth_set_SNPs,
                        truth_set_SNPs_index = truth_set_SNPs_index,
                        training_set_INDELs = genotypeFile.hard_filtered_INDELs_sites_only_vcf,
                        training_set_INDELs_index = genotypeFile.hard_filtered_INDELs_sites_only_vcf_index,
                        training_set_SNPs = genotypeFile.hard_filtered_SNPs_sites_only_vcf,
                        training_set_SNPs_index = genotypeFile.hard_filtered_SNPs_sites_only_vcf_index,
                        container = container_gatk,
                    }
            }

            ## If only one is available, perform partial recalibration combined with hard filtering
            if (!defined(truth_set_SNPs) || !defined(truth_set_INDELs)) {
                call partialVariantRecalibration.partialVQSR as partialVQSR {
                    input:
                        ref = ref,
                        ref_dict = ref_dict,
                        ref_fai = ref_fai,
                        groupName = genotypeFile.groupName,
                        unfiltered_vcf = genotypeFile.unfiltered_vcf,
                        unfiltered_vcf_index = genotypeFile.unfiltered_vcf_index,
                        unfiltered_vcf_sites_only = genotypeFile.unfiltered_sites_only_vcf,
                        unfiltered_vcf_sites_only_index = genotypeFile.unfiltered_sites_only_vcf_index,
                        type = if defined(truth_set_SNPs) && !defined(truth_set_INDELs) then "SNP" else "INDEL",
                        truth_set = if defined(truth_set_SNPs) && !defined(truth_set_INDELs) then truth_set_SNPs else truth_set_INDELs,
                        truth_set_index = if defined(truth_set_SNPs) && !defined(truth_set_INDELs) then truth_set_SNPs_index else truth_set_INDELs_index,
                        training_set = if defined(truth_set_SNPs) && !defined(truth_set_INDELs) then genotypeFile.hard_filtered_SNPs_sites_only_vcf else genotypeFile.hard_filtered_INDELs_sites_only_vcf,
                        training_set_index = if defined(truth_set_SNPs) && !defined(truth_set_INDELs) then genotypeFile.hard_filtered_SNPs_sites_only_vcf_index else genotypeFile.hard_filtered_INDELs_sites_only_vcf_index,
                        container = container_gatk,
                    }
            }

        } ## Stop scattering over genotype sets per group
        
        ## If only one is available, merge the partially VQSR-filtered variants from all groups with the opposite hard-filtered variants from all groups
        if (!defined(truth_set_SNPs) || !defined(truth_set_INDELs)) {

            ## VQSR'd INDELs, hard-filtered SNPs
            if (!defined(truth_set_SNPs) && defined(truth_set_INDELs)) {
                call processVCFs.mergeManyVCFs as mergeINDELsFromVQSRwithHardFilteredSNPs {
                    input:
                        input_array_1_vcf = select_all(partialVQSR.output_vcf),
                        input_array_1_vcf_index = select_all(partialVQSR.output_vcf_index),
                        input_array_2_vcf = select_all(pipe.hard_filtered_SNPs_sites_only_vcf),
                        input_array_2_vcf_index = select_all(pipe.hard_filtered_SNPs_sites_only_vcf_index),
                        ref_dict = ref_dict,
                        output_vcf_title = "mergeINDELsFromVQSRwithHardFilteredSNPs",
                        # Runtime
                        container = container_gatk,
                }              
            }
            ## VQSR'd SNPs, hard-filtered INDELs
            if (defined(truth_set_SNPs) && !defined(truth_set_INDELs)) {
                call processVCFs.mergeManyVCFs as mergeSNPsFromVQSRwithHardFilteredINDELs {
                    input:
                        input_array_1_vcf = select_all(partialVQSR.output_vcf),
                        input_array_1_vcf_index = select_all(partialVQSR.output_vcf_index),
                        input_array_2_vcf = select_all(pipe.hard_filtered_INDELs_sites_only_vcf),
                        input_array_2_vcf_index = select_all(pipe.hard_filtered_INDELs_sites_only_vcf_index),
                        ref_dict = ref_dict,
                        output_vcf_title = "mergeSNPsFromVQSRwithHardFilteredINDELs",
                        # Runtime
                        container = container_gatk,
                }
            }

        }

        ## If both are available, merge the fully VQSR-filtered variants from all groups
        if (defined(truth_set_SNPs) && defined(truth_set_INDELs)) {

            # Merge variants into master loci
            call processVCFs.mergeVCFs as mergeAllFromVQSR {
                input:
                    input_vcf_array = select_all(VQSR.output_vcf),
                    input_vcf_index_array = select_all(VQSR.output_vcf_index),
                    ref_dict = ref_dict,
                    output_vcf_title = "mergeAllFromVQSR",
                    # Runtime
                    container = container_gatk,
            }

        }

    }

    ####################################################################
    ## IF MODE = FINAL: Perform final round of variant calling
    ####################################################################

    if (mode_is_final) {

        File masterLociVcf = select_first([mergeMasterLocifromHF.output_merged_vcf, mergeAllFromVQSR.output_merged_vcf, mergeINDELsFromVQSRwithHardFilteredSNPs.output_merged_vcf, mergeSNPsFromVQSRwithHardFilteredINDELs.output_merged_vcf])
        File masterLociVcfIndex = select_first([mergeMasterLocifromHF.output_merged_vcf_index, mergeAllFromVQSR.output_merged_vcf_index, mergeINDELsFromVQSRwithHardFilteredSNPs.output_merged_vcf_index, mergeSNPsFromVQSRwithHardFilteredINDELs.output_merged_vcf_index])

        ## Split masterLociVCF into interval lists for scattering

        call processVCFs.splitIntervals as splitIntervals {
            input:
                ref = ref,
                ref_dict = ref_dict,
                ref_fai = ref_fai,
                masterLociVcf = masterLociVcf,
                masterLociVcfIndex = masterLociVcfIndex,
                # Runtime
                container = container_gatk,
        }


        ## Scatter over each user-defined scatter interval
        scatter (intervalList in splitIntervals.interval_lists) {

            ## Scatter over each BAM file
            scatter (bamFile in bamList) {

                ## Call haplotypes (i.e. per BAM, per scatter)
                call bamToVcf.haplotypeCaller as haplotypeCallerFinal {
                    input:
                        ref = ref,
                        ref_dict = ref_dict,
                        ref_fai = ref_fai,
                        input_bam = bamFile.bam,
                        input_bam_index = bamFile.bam_index,
                        sampleGroup = "FinalCallset",
                        title_gvcf = bamFile.sampleName + "_haplotypes_final_" + basename(intervalList),
                        masterLociVcf = intervalList,
                        # Runtime
                        container = container_gatk,
                }

            }

            ## Import haplotypes to GenomicsDBs
            call bamToVcf.genomicsDBImport as genomicsDBImportFinal {
                input:
                    input_gvcfs = haplotypeCallerFinal.output_gvcf,
                    input_gvcfs_indexes = haplotypeCallerFinal.output_gvcf_index,
                    title_gdb = "gdb_final_" + basename(intervalList),
                    merge_contigs_into_num_partitions = merge_contigs_into_num_partitions,
                    masterLociVcf = intervalList,
                    # Runtime
                    container = container_gatk,
            }

            ## Genotype the final callset
            call bamToVcf.genotypeGenomicsDB as genotypeGenomicsDBFinal {
                input:
                    ref = ref,
                    ref_dict = ref_dict,
                    ref_fai = ref_fai,
                    input_genomicsdb = genomicsDBImportFinal.output_genomicsdb,
                    input_genomicsdb_title = genomicsDBImportFinal.output_genomicsdb_title,
                    output_vcf_title = "final_variant_calls_" + basename(intervalList) + "_genotyped",
                    groupName = "FinalCallset",
                    masterLociVcf = intervalList,
                    # Runtime
                    container = container_gatk,
            }

        }

        ## Gather the final callsets across scatters/interval lists
        call processVCFs.gatherVCFs as gatherFinalCallset {
            input:
                input_vcfs = genotypeGenomicsDBFinal.output_genotypes_vcf,
                input_vcfs_indexes = genotypeGenomicsDBFinal.output_genotypes_vcf_index,
                output_vcf_title = "FinalCallset",
                # Runtime
                container = container_gatk,
        }

    }

    output {
        ## Outputs from initial and repeat modes
        Array[File]? table_before = BQSR.table_before
        Array[File]? table_after = BQSR.table_after
        Array[File]? plots = BQSR.plots
        Array[File]? recalibrated_bam = BQSR.recalibrated_bam
        Array[File]? recalibrated_bam_index = BQSR.recalibrated_bam_index
        ## Outputs from final mode
        File? masterLoci = masterLociVcf
        File? masterLoci_index = masterLociVcfIndex
        File? finalCallset = gatherFinalCallset.output_vcf
        File? finalCallsetIndex = gatherFinalCallset.output_vcf_index
    }
    
}
