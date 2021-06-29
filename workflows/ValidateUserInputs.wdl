version development

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

    String pipeline_version = "pre-alpha"

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
        Boolean mode_is_initial
        Boolean validate_reference_vcf = false # options: true / false
        Boolean validate_truth_sets = false # options: true / false
        String container_gatk
        String container_python
    }

    ##########################################################################
    ## VALIDATE REFERENCE VCF
    ## The default is false. This can be enabled by setting the Boolean
    ## 'validate_reference_vcf' to true. It will only run in 'initial' mode.
    ##########################################################################

    if (mode_is_initial && validate_reference_vcf) {
        call QC.validateVCF as ValidateReferenceVcf {
            input:
                variant_file = ref,
                ref_dict = ref_dict,
                ref_fai = ref_fai,
                optional_arguments = "-Xtype ALL",
                # Runtime options
                container = container_gatk,
        }
    }

    ##########################################################################
    ## VALIDATE TRUTH SETS IF AVAILABLE
    ## This will always run in initial mode. To force validation in subsequent
    ## modes, set the Boolean 'validate_truth_sets' to true.
    ##########################################################################

    if (mode_is_initial || validate_truth_sets) {
        if (defined(truth_set_INDELs)) {
            call QC.validateVCF as ValidateTruthINDELs {
                input:
                    variant_file = truth_set_INDELs,
                    variant_file_index = truth_set_INDELs_index,
                    ref = ref,
                    ref_dict = ref_dict,
                    ref_fai = ref_fai,
                    optional_arguments = "-Xtype CHR_COUNTS -Xtype IDS -Xtype ALLELES",
                    # Runtime options
                    container = container_gatk,
            }
        }
        if (defined(truth_set_SNPs)) {
            call QC.validateVCF as ValidateTruthSNPs {
                input:
                    variant_file = truth_set_SNPs,
                    variant_file_index = truth_set_SNPs_index,
                    ref = ref,
                    ref_dict = ref_dict,
                    ref_fai = ref_fai,
                    optional_arguments = "-Xtype CHR_COUNTS -Xtype IDS -Xtype ALLELES",
                    # Runtime options
                    container = container_gatk,
            }
        }
    }

    ##########################################################################
    ## VALIDATE THE USER-PROVIDED DATA FOR EACH SAMPLE RECORD
    ## This is performed within this WDL and the ExitWithError task is
    ## called if needed. The rationale for this is to avoid scattering
    ## each sample to a separate task, which would require resource requests
    ## per task and thus be more computationally expensive.
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
                container = container_python,
        }
    }

    ##########################################################################
    ## VALIDATE THE USER-PROVIDED DATA ACROSS THE COHORT
    ## This checks for common errors in the sample and scatter lists.
    ##########################################################################

    call QC.validateCohort as validateCohort {
        input:
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
        Array[String] groupNames = select_all(validateRecords.groupNames)
        Array[String] sampleNames = select_all(validateRecords.sampleNames)
        Array[String] bams = select_all(validateRecords.bams)
        Array[String] bam_indexes = select_all(validateRecords.bam_indexes)
        Array[String] unmapped_bams = select_all(validateRecords.unmapped_bams)
        Boolean multiple_taxonomic_groups = validateCohort.multiple_taxonomic_groups
    }

}
