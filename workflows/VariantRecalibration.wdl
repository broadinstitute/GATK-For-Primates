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
import "../tasks/ProcessVCFs.wdl" as processVCFs
import "../tasks/Recalibration.wdl" as VQSR

##########################################################################
## WORKFLOW DEFINITION
##########################################################################

workflow VQSR {

    input {
        File ref
        File ref_dict
        File ref_fai
        String groupName
        File unfiltered_vcf
        File unfiltered_vcf_index
        File unfiltered_vcf_sites_only
        File unfiltered_vcf_sites_only_index
        File? truth_set_INDELs
        File? truth_set_INDELs_index
        File? truth_set_SNPs
        File? truth_set_SNPs_index
        File training_set_INDELs
        File training_set_INDELs_index
        File training_set_SNPs
        File training_set_SNPs_index
        String container
    }

    ## Recalibrate INDELs
    call VQSR.variantRecalibrator as recalibrateINDELs {
        input:
            ref = ref,
            ref_dict = ref_dict,
            ref_fai = ref_fai,
            groupName = groupName,
            input_vcf_sites_only = unfiltered_vcf_sites_only,
            input_vcf_sites_only_index = unfiltered_vcf_sites_only_index,
            type = "INDEL",
            truth_set = truth_set_INDELs,
            training_set = training_set_INDELs,
            training_set_index = training_set_INDELs_index,
            # Runtime
            container = container,
    }

    ## Recalibrate SNPs
    call VQSR.variantRecalibrator as recalibrateSNPs {
        input:
            ref = ref,
            ref_dict = ref_dict,
            ref_fai = ref_fai,
            groupName = groupName,
            input_vcf_sites_only = unfiltered_vcf_sites_only,
            input_vcf_sites_only_index = unfiltered_vcf_sites_only_index,
            type = "SNP",
            truth_set = truth_set_SNPs,
            training_set = training_set_SNPs,
            training_set_index = training_set_SNPs_index,
            # Runtime
            container = container,
    }

    ## applyVQSR INDELs
    call VQSR.applyVQSR as applyVQSRtoINDELs {
        input:
            ref = ref,
            ref_dict = ref_dict,
            ref_fai = ref_fai,
            groupName = groupName,
            input_vcf = unfiltered_vcf,
            input_vcf_index = unfiltered_vcf_index,
            recalibrated_recal = recalibrateINDELs.output_recalibrated_recal,
            recalibrated_tranches = recalibrateINDELs.output_recalibrated_tranches,
            recalibrated_rscript = recalibrateINDELs.output_recalibrated_rscript,
            exclude_filtered = "false",
            type = "INDEL",
            # Runtime
            container = container,
    }

    ## applyVQSR SNPs
    call VQSR.applyVQSR as applyVQSRtoSNPs {
        input:
            ref = ref,
            ref_dict = ref_dict,
            ref_fai = ref_fai,
            groupName = groupName,
            input_vcf = applyVQSRtoINDELs.output_recalibrated_vcf,
            input_vcf_index = applyVQSRtoINDELs.output_recalibrated_vcf,
            recalibrated_recal = recalibrateSNPs.output_recalibrated_recal,
            recalibrated_tranches = recalibrateSNPs.output_recalibrated_tranches,
            recalibrated_rscript = recalibrateSNPs.output_recalibrated_rscript,
            type = "SNP",
            exclude_filtered = "true",
            # Runtime
            container = container,
    }

   ## Make sites only file -- this is the final callset per group
    call processVCFs.makeSitesOnly as makeSitesOnlyRecals {
        input:
            input_vcf = applyVQSRtoSNPs.output_recalibrated_vcf,
            input_vcf_index = applyVQSRtoSNPs.output_recalibrated_vcf_index,
            groupName = groupName,
            # Runtime
            container = container,
    }

    output {
        File output_vcf = makeSitesOnlyRecals.output_vcf
        File output_vcf_index = makeSitesOnlyRecals.output_vcf_index
    }

}
