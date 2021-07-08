version 1.0

## Copyright Broad Institute and Wisconsin National Primate Research Center,
## University of Wisconsin-Madison, 2021
## 
## Tasks from the complete germline short variant discovery pipeline
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

##########################################################################
## *** TASK: haplotypeCaller ***
##########################################################################
## Calls haplotypes from input BAM files. Optionally adds padding
## if this is the final call to re-genotype at the master list of loci.
##########################################################################

task haplotypeCaller {
    input {
        File ref
        File ref_dict
        File ref_fai
        File input_bam
        File input_bam_index
        String title_gvcf
        String? scatterIntervals
        File? masterLociVcf
        File? masterLociVcfIndex
        String sampleGroup
        # runtime
        String container
        Int? runtime_set_preemptible_tries
        Int? runtime_set_cpu
        Int? runtime_set_memory
        Int? runtime_set_disk
        Int? runtime_set_max_retries
        Boolean use_ssd = false

    }
    ## Runtime parameters
    Float size_input_files = size(ref, "GB") + size(ref_dict, "GB")  + size(ref_fai, "GB") + size(input_bam, "GB") + size(input_bam_index, "GB") + size(masterLociVcf, "GB") + size(masterLociVcfIndex, "GB")
    Int runtime_calculated_disk = ceil(size_input_files * 2.5)
    Int command_mem_gb = select_first([runtime_set_memory, 7]) - 1
    ## Task-specific parameters
    String? add_padding = if defined(masterLociVcf) then "--interval-padding 100" else ""
    command {
        gatk \
        HaplotypeCaller --java-options "-Xmx~{command_mem_gb}G" \
        -R ~{ref} \
        -I ~{input_bam} \
        -O ~{title_gvcf}.g.vcf.gz \
        -L ~{scatterIntervals}~{masterLociVcf} \
        -ERC GVCF \
        -G AS_StandardAnnotation \
        -GQB 10 -GQB 20 -GQB 30 -GQB 40 -GQB 50 -GQB 60 -GQB 70 -GQB 80 -GQB 90 \
        ~{add_padding}
    }
    output {
        File output_gvcf = "~{title_gvcf}" + ".g.vcf.gz"
        File output_gvcf_index = "~{title_gvcf}" + ".g.vcf.gz.tbi"
        String output_sampleGroup = "~{sampleGroup}"
        String output_title_gvcf = "~{title_gvcf}"
    }
    runtime {
        docker: container
        cpu: select_first([runtime_set_cpu, 1])
        gpu: false
        memory: select_first([runtime_set_memory, 7]) + " GB"
        disks: "local-disk " + select_first([runtime_set_disk, runtime_calculated_disk]) + if use_ssd then " SSD" else " HDD"
        maxRetries: select_first([runtime_set_max_retries, 0])
        preemptible: select_first([runtime_set_preemptible_tries, 5])
        returnCodes: 0
     }
}


##########################################################################
## *** TASK: genomicsDBImport ***
##########################################################################
## Imports haplotypes into genomicsDBs.
##########################################################################

task genomicsDBImport {
    input {
        Array[File] input_gvcfs
        Array[File] input_gvcfs_indexes
        String title_gdb
        String? scatterIntervals
        File? masterLociVcf
        File? masterLociVcfIndex
        Int? merge_contigs_into_num_partitions
        # runtime
        String container
        Int? runtime_set_preemptible_tries
        Int? runtime_set_cpu
        Int? runtime_set_memory
        Int? runtime_set_disk
        Int? runtime_set_max_retries
        Boolean use_ssd = false
    }
    ## Runtime parameters
    # This requires attention
    Float size_input_files = size(input_gvcfs, "GB") + size(input_gvcfs_indexes, "GB") ##+ size(masterLociVcf, "GB") + size(masterLociVcfIndex, "GB")
    Int runtime_calculated_disk = ceil(size_input_files * 2.5)
    Int command_mem_gb = select_first([runtime_set_memory, 9])
    ## Task-specific parameters
    ## This requires attention:
    ## Int merge_contigs_value = if defined(merge_contigs_into_num_partitions) then merge_contigs_into_num_partitions else "0"
    String? add_padding = if defined(masterLociVcf) then "--interval-padding 100" else ""
    command <<<
    set -euo pipefail

        echo "~{title_gdb}"
        echo "X"

        gatk \
        GenomicsDBImport --java-options "-Xmx~{command_mem_gb}G" \
        -V ~{sep=' -V ' input_gvcfs} \
        --genomicsdb-workspace-path ~{title_gdb} \
        --batch-size 50 \
        --reader-threads 5 \
        --consolidate \
        -L ~{scatterIntervals}~{masterLociVcf} \
        ~{add_padding}

        tar -cf ~{title_gdb}.tar ~{title_gdb}
    >>>
    output {
        File output_genomicsdb = "~{title_gdb}.tar"
        String output_genomicsdb_title = "~{title_gdb}"
    }
    runtime {
        docker: container
        cpu: select_first([runtime_set_cpu, 1])
        gpu: false
        memory: select_first([runtime_set_memory, 8]) + " GB"
        disks: "local-disk " + select_first([runtime_set_disk, runtime_calculated_disk, 8]) + if use_ssd then " SSD" else " HDD"
        maxRetries: select_first([runtime_set_max_retries, 0])
        preemptible: select_first([runtime_set_preemptible_tries, 1])
        returnCodes: 0
     }
}


##########################################################################
## *** TASK: genotypeGenomicsDB ***
##########################################################################
## Genotypes input genomicsDBs. Optionally calls only at master list of
## loci, including non-variant sites, if this is the final call.
##########################################################################

task genotypeGenomicsDB {
    input {
        File ref
        File ref_dict
        File ref_fai
        String? groupName
        File input_genomicsdb
        String input_genomicsdb_title
        String output_vcf_title
        File? masterLociVcf
        File? masterLociVcfIndex
        String? scatterName        
        String? scatterIntervals
        # runtime
        String container
        Int? runtime_set_preemptible_tries
        Int? runtime_set_cpu
        Int? runtime_set_memory
        Int? runtime_set_disk
        Int? runtime_set_max_retries
        Boolean use_ssd = false
    }
    ## Runtime parameters
    Float size_input_files = size(ref, "GB") + size(ref_dict, "GB") + size(ref_fai, "GB") + size(input_genomicsdb, "GB") ##+ size(masterLociVcf, "GB") + size(masterLociVcfIndex, "GB")
    Int runtime_calculated_disk = ceil(size_input_files * 2.5) + 2
    Int command_mem_gb = select_first([runtime_set_memory, 8]) - 1
    ## Task-specific parameters
    String? genotype_non_variant_sites = if defined(masterLociVcf) then "--include-non-variant-sites" else ""
    command <<<
    set -euo pipefail

        tar -xf ~{input_genomicsdb}

        gatk \
        GenotypeGVCFs --java-options "-Xmx~{command_mem_gb}G" \
        -R ~{ref} \
        -V gendb://~{input_genomicsdb_title} \
        -L ~{scatterIntervals}~{masterLociVcf} \
        --only-output-calls-starting-in-intervals \
        -O ~{output_vcf_title}.vcf.gz \
        -G StandardAnnotation -G AS_StandardAnnotation \
        ~{genotype_non_variant_sites}
        
    >>>
    output {
        File output_genotypes_vcf = "~{output_vcf_title}.vcf.gz"
        File output_genotypes_vcf_index = "~{output_vcf_title}.vcf.gz.tbi"
        String output_groupName = "~{groupName}"
        String output_scatterName = "~{scatterName}"
    }
    runtime {
        docker: container
        cpu: select_first([runtime_set_cpu, 1])
        gpu: false
        memory: select_first([runtime_set_memory, 8]) + " GB"
        disks: "local-disk " + select_first([runtime_set_disk, runtime_calculated_disk]) + if use_ssd then " SSD" else " HDD"
        maxRetries: select_first([runtime_set_max_retries, 0])
        preemptible: select_first([runtime_set_preemptible_tries, 5])
        returnCodes: 0
     }
}
