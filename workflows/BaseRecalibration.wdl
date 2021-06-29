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
import "../tasks/Recalibration.wdl" as BQSR

##########################################################################
## WORKFLOW DEFINITION
##########################################################################

workflow BQSR {

    input {
        File ref
        File ref_dict
        File ref_fai
        File input_bam
        File input_bam_index
        Array[File] input_SNP_sites
        Array[File] input_SNP_sites_indexes
        Array[File] input_INDEL_sites
        Array[File] input_INDEL_sites_indexes
        String sampleName
        String container_gatk
        String container_gitc
        String path_to_gitc_gatk
    }

    ## Call baseRecalibrator etc. for each BAM file
    call BQSR.baseRecalibrator as baseRecalibrator {
    input:
            ref = ref,
            ref_dict = ref_dict,
            ref_fai = ref_fai,
            input_bam = input_bam,
            input_bam_index = input_bam_index,
            input_SNP_sites = input_SNP_sites,
            input_SNP_sites_indexes = input_SNP_sites_indexes,
            input_INDEL_sites = input_INDEL_sites,
            input_INDEL_sites_indexes = input_INDEL_sites_indexes,
            sampleName = sampleName,
            # Runtime
            container = container_gatk,

    }
    
    ## Make plots for each BAM file -- requires separate task as R and ggplot needed
    ## and the genomesinthecloud docker uses an old version of GATK
    call BQSR.analyzeCovariates as analyzeCovariates {
    input:
            sampleName = sampleName,
            table_before = baseRecalibrator.table_before,
            table_after = baseRecalibrator.table_after,
            # Runtime
            container = container_gitc,
            path_to_gitc_gatk = path_to_gitc_gatk,
    }

    output {
        File table_before = baseRecalibrator.table_before
        File table_after = baseRecalibrator.table_after
        File recalibrated_bam = baseRecalibrator.recalibrated_bam
        File recalibrated_bam_index = baseRecalibrator.recalibrated_bam_index
        File plots = analyzeCovariates.plots
    }

}
