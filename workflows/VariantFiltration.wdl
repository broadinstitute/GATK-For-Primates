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

##########################################################################
## WORKFLOW DEFINITION
##########################################################################

workflow variantFiltration {

    input {
        File ref
        File ref_dict
        File ref_fai
        File input_vcf
        File input_vcf_index
        String container
    }

	## Split multi-allelics, variant types, hard filter and make sites only
    call processVCFs.splitVariants as splitVariants {
        input:
            ref = ref,
            ref_dict = ref_dict,
            ref_fai = ref_fai,
            input_vcf = input_vcf,
            input_vcf_index = input_vcf_index,
            # Runtime
            container = container,
    }

    ## Hard filter SNPs and make sites only
    call processVCFs.hardFilter as hardFilterSNPs {
        input:
        	ref = ref,
        	ref_dict = ref_dict,
        	ref_fai = ref_fai,
            makeSitesOnly = true,
            input_vcf = splitVariants.output_SNPs,
            input_vcf_index = splitVariants.output_SNPs_index,
            type = "SNPs",
            filters = "-filter 'QD < 2.0' --filter-name 'QD2' -filter 'QUAL < 30.0' --filter-name 'QUAL30' -filter 'SOR > 3.0' --filter-name 'SOR3' -filter 'FS > 60.0' --filter-name 'FS60' -filter 'MQ < 40.0' --filter-name 'MQ40' -filter 'MQRankSum < -12.5' --filter-name 'MQRankSum-12.5' -filter 'ReadPosRankSum < -8.0' --filter-name 'ReadPosRankSum-8'",
            # Runtime
            container = container,
    }

    ## Hard filter INDELs and make sites only
    call processVCFs.hardFilter as hardFilterINDELs {
        input:
        	ref = ref,
        	ref_dict = ref_dict,
        	ref_fai = ref_fai,
            makeSitesOnly = true,
            input_vcf = splitVariants.output_INDELs,
            input_vcf_index = splitVariants.output_INDELs_index,
			type = "INDELs",
            filters = "-filter 'QD < 2.0' --filter-name 'QD2' -filter 'QUAL < 30.0' --filter-name 'QUAL30' -filter 'FS > 200.0' --filter-name 'FS200' -filter 'ReadPosRankSum < -20.0' --filter-name 'ReadPosRankSum-20'",
            # Runtime
            container = container,
    }


    output {
        File SNP_sites_vcf = hardFilterSNPs.output_sites_only_vcf
        File SNP_sites_vcf_index = hardFilterSNPs.output_sites_only_vcf_index
        File INDEL_sites_vcf = hardFilterINDELs.output_sites_only_vcf
        File INDEL_sites_vcf_index = hardFilterINDELs.output_sites_only_vcf_index
    }

}
