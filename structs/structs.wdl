version development

## Copyright Broad Institute and Wisconsin National Primate Research Center,
## University of Wisconsin-Madison, 2021
## 
## Structs from the complete germline short variant discovery pipeline
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

struct sampleInfo {
    String name
    String? RG_ID
    String? RG_LB
    String? RG_SM
    String? RG_PU
    String taxon_group
    File? R1
    File? R2
    File? unmapped_bam
    File? bam
    File? bam_index
}

struct scatterInfo {
    String name
    String intervals
}

struct bamInfo {
    String sampleGroup
    String sampleName
    String bam
    String bam_index
}

struct haplotypeInfo {
    String groupName
    String title_gvcf
    String gvcf
    String gvcf_index
}

struct genotypeInfo {
    String groupName
    String unfiltered_vcf
    String unfiltered_vcf_index
    String unfiltered_sites_only_vcf
    String unfiltered_sites_only_vcf_index
    String hard_filtered_SNPs_vcf
    String hard_filtered_SNPs_vcf_index
    String hard_filtered_INDELs_vcf
    String hard_filtered_INDELs_vcf_index
    String hard_filtered_SNPs_sites_only_vcf
    String hard_filtered_SNPs_sites_only_vcf_index
    String hard_filtered_INDELs_sites_only_vcf
    String hard_filtered_INDELs_sites_only_vcf_index
}
