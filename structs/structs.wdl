version 1.0

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

struct readInfo {
    String read_id
    String sample_id
    String taxon_group
    String? R1
    String? R2
    String? unmapped_bam
    String RG_ID
    String RG_LB
    String RG_SM
    String RG_PU
    String? RG_CN
    String? RG_DS
    String? RG_DT
    String? RG_FO
    String? RG_KS
    String? RG_PG
    String? RG_PI
    String? RG_PM

struct sampleInfo {
    String sample_id
    String taxon_group
    File? bam
    File? bam_index
    File? recalibrated_bam
    File? recalibrated_bam_index
    
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
    String gvcf
}

struct polymorphicRegionsInfo {
    String scatterName
    String intervalList
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
