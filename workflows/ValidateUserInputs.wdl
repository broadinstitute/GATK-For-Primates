version 1.0

## Copyright Broad Institute and Wisconsin National Primate Research Center,
## University of Wisconsin-Madison, 2021
## 
## Sub-workflow from the complete germline short variant discovery pipeline
## optimized for non-human primates. For requirements, expectations and
## outputs, please review the complete documentation at:
## https://github.com/broadinstitute/GATK-For-Primates
##
## LICENSING : 
## This script is released under the WDL source code license (BSD-3) (see LICENSE in 
## https://github.com/broadinstitute/wdl). However, the programs it calls may be
## subject to different licenses. Users are responsible for checking that they are
## authorized to run all programs before running this script. Please see the docker
## containers for detailed licensing information pertaining to the included programs.

import "../structs/structs.wdl"
import "../tasks/QualityControl.wdl" as QC
import "../tasks/Utilities.wdl" as utilities

##########################################################################
## WORKFLOW DEFINITION
##########################################################################

workflow validateUserInputs {

    ##########################################################################
    ## WORKFLOW INPUTS
    ##########################################################################

    input {
        File ref
        File ref_dict
        File ref_fai
        Array[sampleInfo]+ sampleList
        Array[scatterInfo]+ scatterList
        File? truth_set_SNPs
        File? truth_set_SNPs_index
        File? truth_set_INDELs
        File? truth_set_INDELs_index
        File? packaged_polymorphic_regions
        String mode
        Boolean validate_truth_sets = true # options: true / false
        String container_gatk
        String container_python
    }

    ##########################################################################
    ## Fail workflow if `mode` is not 'initial', 'repeat' or 'final'
    ##########################################################################

    if ((mode != "initial") && (mode != "repeat") && (mode != "final")) {
        call QC.failWithError as failWrongMode {
            input:
                message = "Input mode must be either 'initial', 'repeat' or 'final'.",
                # Runtime options
                container = container_python,
        }
    }

    ##########################################################################
    ## Fail workflow in 'repeat'/'final' modes if required input file
    ## `packaged_polymorphic_regions` is missing
    ##########################################################################

    if ((mode != "initial") && (!defined(packaged_polymorphic_regions))) {
        call QC.failWithError as failMissingPackage {
            input:
                message = "Packaged polymorphic regions from 'initial' mode must be provided in 'repeat' or 'final' modes.",
                # Runtime options
                container = container_python,
        }
    }

    ##########################################################################
    ## Unpackage polymorphic region JSON and interval list files, or make blank
    ## JSON in 'initial' mode (as Cromwell dislikes optional inputs/outputs)
    ##########################################################################

    call utilities.unpackagePolymorphicRegions as unpackagePolymorphicRegions {
        input:
            mode = mode,
            packaged_polymorphic_regions = packaged_polymorphic_regions,
            # Runtime options
            container = container_python,
    }

    ##########################################################################
    ## Validate truth sets in 'final' mode only.
    ## This is 'true' by default in the main workflow options
    ##########################################################################

    if ((mode == "final") && (validate_truth_sets)) {
        if (defined(truth_set_INDELs)) {
            call QC.validateVCF as ValidateTruthINDELs {
                input:
                    ref = ref,
                    ref_dict = ref_dict,
                    ref_fai = ref_fai,
                    variant_file = truth_set_INDELs,
                    variant_file_index = truth_set_INDELs_index,
                    optional_arguments = "-Xtype CHR_COUNTS -Xtype IDS -Xtype ALLELES",
                    # Runtime options
                    container = container_gatk,
            }
        }
        if (defined(truth_set_SNPs)) {
            call QC.validateVCF as ValidateTruthSNPs {
                input:
                    ref = ref,
                    ref_dict = ref_dict,
                    ref_fai = ref_fai,
                    variant_file = truth_set_SNPs,
                    variant_file_index = truth_set_SNPs_index,
                    optional_arguments = "-Xtype CHR_COUNTS -Xtype IDS -Xtype ALLELES",
                    # Runtime options
                    container = container_gatk,
            }
        }
    }

    ##########################################################################
    ## Validate user-provided data for each sample record
    ## Validation prevents computationally expensive failures later on!
    ##########################################################################

    scatter (sample in sampleList) {
        call QC.validateRecords as validateRecords {
            input:
                groupName = if defined(sample.taxon_group) then sample.taxon_group else "NULL",
                sampleName = if defined(sample.name) then sample.name else "NULL",
                R1 = if defined(sample.R1) then sample.R1 else "NULL",
                R2 = if defined(sample.R2) then sample.R2 else "NULL",
                RG_ID = if defined(sample.RG_ID) then sample.RG_ID else "NULL",
                RG_SM = if defined(sample.RG_SM) then sample.RG_SM else "NULL",
                RG_LB = if defined(sample.RG_LB) then sample.RG_LB else "NULL",
                RG_PU = if defined(sample.RG_PU) then sample.RG_PU else "NULL",
                bam = if defined(sample.bam) then sample.bam else "NULL",
                bam_index = if defined(sample.bam_index) then sample.bam_index else "NULL",
                unmapped_bam = if defined(sample.unmapped_bam) then sample.unmapped_bam else "NULL",
                # Runtime options
                container = container_python,
        }
    }

    ##########################################################################
    ## Validate user-provided data across the cohort
    ## This checks for common errors in the sample and scatter lists
    ##########################################################################

    scatter (scttr in scatterList) {
        String loop_scatterNames = scttr.name
        String loop_scatterIntervals = scttr.intervals
    }

    call QC.validateCohort as validateCohort {
        input:
            scatterNames = loop_scatterNames,
            scatterIntervals = loop_scatterIntervals,
            groupNames = validateRecords.groupNames,
            sampleNames = validateRecords.sampleNames,
            RG_IDs = validateRecords.RG_IDs,
            RG_SMs = validateRecords.RG_SMs,
            RG_LBs = validateRecords.RG_LBs,
            container = container_python,
    }


    ##########################################################################
    ## WORKFLOW OUTPUTS
    ##########################################################################

    output {
        File polymorphic_regions_json = unpackagePolymorphicRegions.polymorphicRegionsJSON
        Array[File] polymorphicRegionsIntervalLists = unpackagePolymorphicRegions.polymorphicRegionsIntervalLists
        Array[String] groupNames = select_all(validateRecords.groupNames)
        Array[String] sampleNames = select_all(validateRecords.sampleNames)
        Array[String] bams = select_all(validateRecords.bams)
        Array[String] bam_indexes = select_all(validateRecords.bam_indexes)
        Array[String] unmapped_bams = select_all(validateRecords.unmapped_bams)
        Array[String] scatterNames = select_all(validateCohort.output_scatterNames)
        Boolean multiple_taxonomic_groups = validateCohort.multiple_taxonomic_groups
    }

}
