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
## *** TASK: makeSitesOnly ***
##########################################################################
## Makes a sites-only version of the input VCF.
##########################################################################

task makeSitesOnly {
    input {
        File input_vcf
        File input_vcf_index
        String groupName
        # Runtime
        String container
        Int? runtime_set_preemptible_tries
        Int? runtime_set_cpu
        Int? runtime_set_memory
        Int? runtime_set_disk
        Int? runtime_set_max_retries
        Boolean use_ssd = false
    }
    Float size_input_files = size(input_vcf, "GB") + size(input_vcf_index, "GB") 
    Int runtime_calculated_disk = ceil(size_input_files * 2)
    Int command_mem_gb = select_first([runtime_set_memory, 4]) - 1
    String basename = basename(input_vcf, ".vcf.gz")
    command <<<
    
        gatk MakeSitesOnlyVcf --java-options "-Xmx~{command_mem_gb}G" \
        -I ~{input_vcf} \
        -O ~{basename}_sitesOnly.vcf.gz

    >>>
    output {
        File output_vcf = "~{basename}_sitesOnly.vcf.gz"
        File output_vcf_index = "~{basename}_sitesOnly.vcf.gz.tbi"
        String output_groupName = "~{groupName}"
    }
    runtime {
        docker: container
        cpu: select_first([runtime_set_cpu, 1])
        gpu: false
        memory: select_first([runtime_set_memory, 4]) + " GB"
        disks: "local-disk " + select_first([runtime_set_disk, runtime_calculated_disk]) + if use_ssd then " SSD" else " HDD"
        maxRetries: select_first([runtime_set_max_retries, 0])
        preemptible: select_first([runtime_set_preemptible_tries, 5])
        returnCodes: 0
     }
}

##########################################################################
## *** TASK: mergeVCFs ***
##########################################################################
## Merges multiple VCFs into a single VCF.
##########################################################################

task mergeVCFs {
    input {
        Array[File]? input_vcf_array
        Array[File]? input_vcf_index_array
        File? input_vcf
        File? input_vcf_index
        File ref_dict
        String output_vcf_title
        String? groupName
        # Runtime
        String container
        Int? runtime_set_preemptible_tries
        Int? runtime_set_cpu
        Int? runtime_set_memory
        Int? runtime_set_disk
        Int? runtime_set_max_retries
        Boolean use_ssd = false
    }
    Float size_input_files = size(input_vcf_array, "GB") + size(input_vcf_index_array, "GB") + size(input_vcf, "GB") + size(input_vcf_index, "GB") + size(ref_dict, "GB")
    Int runtime_calculated_disk = ceil(size_input_files * 2.5)
    Int command_mem_gb = select_first([runtime_set_memory, 7]) - 1
    String param_array = if defined(input_vcf_array) then "-I ~{sep=" -I " input_vcf_array}" else ""
    String param_file = if defined(input_vcf) then "-I ~{input_vcf}" else ""
    command <<<

        gatk \
        MergeVcfs --java-options "-Xmx~{command_mem_gb}G" \
        ~{param_array} ~{param_file} \
        -O ~{output_vcf_title}.vcf.gz \
        -D ~{ref_dict}

    >>>
    output {
        File output_merged_vcf = "~{output_vcf_title}.vcf.gz"
        File output_merged_vcf_index = "~{output_vcf_title}.vcf.gz.tbi"
        String output_groupName = "~{groupName}"
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
## *** TASK: mergeManyVCFs ***
##########################################################################
## Merges many VCFs (presented in two Array[File])
##########################################################################

task mergeManyVCFs {
    input {
        File ref_dict
        Array[File] input_array_1_vcf
        Array[File] input_array_1_vcf_index
        Array[File] input_array_2_vcf
        Array[File] input_array_2_vcf_index
        String output_vcf_title
        # Runtime
        String container
        Int? runtime_set_preemptible_tries
        Int? runtime_set_cpu
        Int? runtime_set_memory
        Int? runtime_set_disk
        Int? runtime_set_max_retries
        Boolean use_ssd = false
    }
    Float size_input_files = size(input_array_1_vcf, "GB") + size(input_array_1_vcf_index, "GB") + size(input_array_2_vcf, "GB") + size(input_array_2_vcf_index, "GB") + size(ref_dict, "GB")
    Int runtime_calculated_disk = ceil(size_input_files * 2.5)
    Int command_mem_gb = select_first([runtime_set_memory, 6]) - 1
    command <<<

        gatk \
        MergeVcfs --java-options "-Xmx~{command_mem_gb}G" \
        -I ~{sep=' -I ' input_array_1_vcf} -I ~{sep=' -I ' input_array_2_vcf} \
        -O ~{output_vcf_title}.vcf.gz \
        -D ~{ref_dict}

    >>>
    output {
        File output_merged_vcf = "~{output_vcf_title}.vcf.gz"
        File output_merged_vcf_index = "~{output_vcf_title}.vcf.gz.tbi"
    }
    runtime {
        docker: container
        cpu: select_first([runtime_set_cpu, 1])
        gpu: false
        memory: select_first([runtime_set_memory, 6]) + " GB"
        disks: "local-disk " + select_first([runtime_set_disk, runtime_calculated_disk]) + if use_ssd then " SSD" else " HDD"
        maxRetries: select_first([runtime_set_max_retries, 0])
        preemptible: select_first([runtime_set_preemptible_tries, 5])
        returnCodes: 0
     }
}

##########################################################################
## *** TASK: polymorphicSitesToRegions ***
##########################################################################
## Merges VCFs then pre-processes into interval lists
##########################################################################

task polymorphicSitesToRegions {
    input {
        File ref
        File ref_dict
        File ref_fai
        Array[File] input_vcf_array
        Array[File] input_vcf_index_array
        String scatterName
        # Runtime
        String container
        Int? runtime_set_preemptible_tries
        Int? runtime_set_cpu
        Int? runtime_set_memory
        Int? runtime_set_disk
        Int? runtime_set_max_retries
        Boolean use_ssd = false
    }
    Float size_input_files = size(ref, "GB") + size(ref_dict, "GB") + size(ref_fai, "GB") + size(input_vcf_array, "GB") + size(input_vcf_index_array, "GB")
    Int runtime_calculated_disk = ceil(size_input_files * 2.5)
    Int command_mem_gb = select_first([runtime_set_memory, 7]) - 1
    command <<<

        gatk \
        MergeVcfs --java-options "-Xmx~{command_mem_gb}G" \
        -I ~{sep=" -I " input_vcf_array} \
        -O ~{scatterName}_polymorphic_sites.vcf.gz \
        -D ~{ref_dict}

        gatk \
        VcfToIntervalList --java-options "-Xmx~{command_mem_gb}G" \
        -I ~{scatterName}_polymorphic_sites.vcf.gz \
        -O ~{scatterName}_polymorphic_sites.interval_list \
        
        gatk \
        PreprocessIntervals --java-options "-Xmx~{command_mem_gb}G" \
        -R ~{ref} \
        -L ~{scatterName}_polymorphic_sites.interval_list \
        --bin-length 0 \
        --padding 250 \
        --interval-merging-rule OVERLAPPING_ONLY \
        -O ~{scatterName}_polymorphic_regions.interval_list

    >>>
    output {
        String output_scatterName = "~{scatterName}"
        File output_intervalList = "~{scatterName}_polymorphic_regions.interval_list"
        String output_intervalList_filename = "~{scatterName}_polymorphic_regions.interval_list"
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
## *** TASK: splitVariants ***
##########################################################################
## Splits multi-allelic rows without trimming alleles, then splits again
## into separate SNP and INDEL files.
##########################################################################

task splitVariants {
    input {
        File ref
        File ref_dict
        File ref_fai
        File input_vcf
        File input_vcf_index
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
    Float size_input_files = size(ref, "GB") + size(ref_dict, "GB") + size(ref_fai, "GB") + size(input_vcf, "GB") + size(input_vcf_index, "GB")
    Int runtime_calculated_disk = ceil(size_input_files * 2.5)
    Int command_mem_gb = select_first([runtime_set_memory, 4]) - 1
    String basename = basename("~{input_vcf}", ".vcf.gz")
    command <<<

    gatk \
    LeftAlignAndTrimVariants --java-options "-Xmx~{command_mem_gb}G" \
    -R ~{ref} \
    -V ~{input_vcf} \
    -O ~{basename}_biallelic.vcf.gz \
    --split-multi-allelics \
    --dont-trim-alleles

    gatk \
    SplitVcfs --java-options "-Xmx~{command_mem_gb}G" \
    -I ~{basename}_biallelic.vcf.gz \
    --SNP_OUTPUT ~{basename}_biallelic_SNPs.vcf.gz \
    --INDEL_OUTPUT ~{basename}_biallelic_INDELs.vcf.gz \
    --STRICT false

    >>>
    output {
        File output_SNPs = "~{basename}_biallelic_SNPs.vcf.gz"
        File output_SNPs_index = "~{basename}_biallelic_SNPs.vcf.gz.tbi"
        File output_INDELs = "~{basename}_biallelic_INDELs.vcf.gz"
        File output_INDELs_index = "~{basename}_biallelic_INDELs.vcf.gz.tbi"
    }
    runtime {
        docker: container
        cpu: select_first([runtime_set_cpu, 1])
        gpu: false
        memory: select_first([runtime_set_memory, 4]) + " GB"
        disks: "local-disk " + select_first([runtime_set_disk, runtime_calculated_disk]) + if use_ssd then " SSD" else " HDD"
        maxRetries: select_first([runtime_set_max_retries, 0])
        preemptible: select_first([runtime_set_preemptible_tries, 5])
        returnCodes: 0
     }
}

##########################################################################
## *** TASK: hardFilter ***
##########################################################################
## Hard-filters input VCF file according to provided parameters.
##########################################################################

task hardFilter {
    input {
        File ref
        File ref_dict
        File ref_fai
        File input_vcf
        File input_vcf_index
        String filters
        String? optional_parameters
        String? type
        Boolean makeSitesOnly
        # runtime
        String container
        Int? runtime_set_preemptible_tries
        Int? runtime_set_cpu
        Int? runtime_set_memory
        Int? runtime_set_disk
        Int? runtime_set_max_retries
        Boolean use_ssd = false
    }
    Float size_input_files = size(ref, "GB") + size(ref_dict, "GB") + size(ref_fai, "GB") + size(input_vcf, "GB") + size(input_vcf_index, "GB")
    Int runtime_calculated_disk = ceil(size_input_files * 2.5) + 2
    Int command_mem_gb = select_first([runtime_set_memory, 7]) - 1
    String basename = basename(input_vcf, ".vcf.gz")
    String callMakesSitesOnlyVcf = if makeSitesOnly then "gatk MakeSitesOnlyVcf -I ~{basename}_hardFiltered.vcf.gz -O ~{basename}_hardFiltered_sitesOnly.vcf.gz" else ""
    command <<<

        gatk \
        VariantFiltration --java-options "-Xmx~{command_mem_gb}G" \
        -R ~{ref} \
        -V ~{input_vcf} \
        -O ~{basename}_temp.vcf.gz \
        ~{filters} ~{optional_parameters} \
        --verbosity ERROR

        gatk \
        SelectVariants --java-options "-Xmx~{command_mem_gb}G" \
        -R ~{ref} \
        -V ~{basename}_temp.vcf.gz \
        -O ~{basename}_hardFiltered.vcf.gz \
        --exclude-filtered
        
        touch ~{basename}_hardFiltered_sitesOnly.vcf.gz
        touch ~{basename}_hardFiltered_sitesOnly.vcf.gz.tbi

        ~{callMakesSitesOnlyVcf}

    >>>
    output {
        File output_vcf = "~{basename}_hardFiltered.vcf.gz"
        String output_vcf_string = "~{basename}_hardFiltered.vcf.gz"
        File output_vcf_index = "~{basename}_hardFiltered.vcf.gz.tbi"
        File output_sites_only_vcf = "~{basename}_hardFiltered_sitesOnly.vcf.gz"
        File output_sites_only_vcf_index = "~{basename}_hardFiltered_sitesOnly.vcf.gz.tbi"
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
## *** TASK: gatherVCFs ***
##########################################################################
## Gathers non-overlapping VCF files into a single, sorted VCF file.
##########################################################################

task gatherVCFs {
    input {
        Array[File] input_vcfs
        Array[File] input_vcfs_indexes
        String output_vcf_title
        String? groupName
        # Runtime
        String container
        Int? runtime_set_preemptible_tries
        Int? runtime_set_cpu
        Int? runtime_set_memory
        Int? runtime_set_disk
        Int? runtime_set_max_retries
        Boolean use_ssd = false
    }
    Float size_input_files = size(input_vcfs, "GB") + size(input_vcfs_indexes, "GB")
    Int runtime_calculated_disk = ceil(size_input_files * 2.5)
    Int command_mem_gb = select_first([runtime_set_memory, 7]) - 1
    command <<<
    
            gatk \
            GatherVcfs --java-options "-Xmx~{command_mem_gb}G" \
            -I ~{sep=" -I " input_vcfs} \
            --REORDER_INPUT_BY_FIRST_VARIANT
            -O ~{output_vcf_title}.vcf.gz \


            ## Index creation from block-compressed files not yet supported
            ## See https://github.com/broadinstitute/picard/issues/789
            gatk IndexFeatureFile \
            -I ~{output_vcf_title}.vcf.gz

    >>>
    output {
        File output_vcf = "~{output_vcf_title}.vcf.gz"
        File output_vcf_index = "~{output_vcf_title}.vcf.gz.tbi"
        String output_groupName = "~{groupName}"
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
