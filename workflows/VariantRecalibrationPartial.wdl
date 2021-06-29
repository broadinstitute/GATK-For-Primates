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
import "../tasks/ProcessVCFs.wdl" as processVCFs
import "../tasks/Recalibration.wdl" as VQSR

##########################################################################
## WORKFLOW DEFINITION
##########################################################################

workflow partialVQSR {

    input {
        File ref
        File ref_dict
        File ref_fai
        String groupName
        File unfiltered_vcf
        File unfiltered_vcf_index
        File unfiltered_vcf_sites_only
        File unfiltered_vcf_sites_only_index
        String type
        File? truth_set
        File? truth_set_index
        File training_set
        File training_set_index
        String container
    }

    ## Recalibrate the available one
    call VQSR.variantRecalibrator as recalibratePartial {
        input:
            ref = ref,
            ref_dict = ref_dict,
            ref_fai = ref_fai,
            groupName = groupName,
            input_vcf_sites_only = unfiltered_vcf_sites_only,
            input_vcf_sites_only_index = unfiltered_vcf_sites_only_index,
            type = type,
            truth_set = truth_set,
            truth_set_index = truth_set_index,
            training_set = training_set,
            training_set_index = training_set_index,
            # Runtime
            container = container,
    }

    ## applyVQSR to the available one
    call VQSR.applyVQSR as applyVQSRtoPartial {
        input:
            ref = ref,
            ref_dict = ref_dict,
            ref_fai = ref_fai,
            groupName = groupName,
            input_vcf = unfiltered_vcf,
            input_vcf_index = unfiltered_vcf_index,
            recalibrated_recal = recalibratePartial.output_recalibrated_recal,
            recalibrated_tranches = recalibratePartial.output_recalibrated_tranches,
            recalibrated_rscript = recalibratePartial.output_recalibrated_rscript,
            type = type,
            exclude_filtered = true,
            # Runtime
            container = container,
    }

   ## Make sites only of the recalibrated available set
    call processVCFs.makeSitesOnly as makeSitesOnlyPartials {
        input:
            input_vcf = applyVQSRtoPartial.output_recalibrated_vcf,
            input_vcf_index = applyVQSRtoPartial.output_recalibrated_vcf_index,
            groupName = groupName,
            output_vcf_title = "recalibrated_available_sites_only.",
            # Runtime
            container = container,
    }


    output {
        File output_vcf = makeSitesOnlyPartials.output_vcf
        File output_vcf_index = makeSitesOnlyPartials.output_vcf_index
    }

}
