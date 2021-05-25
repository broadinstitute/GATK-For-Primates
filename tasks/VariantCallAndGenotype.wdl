version development

task HaplotypeCaller {
    input {
        File ref
        File ref_dict
        Array[File]+ ref_idxs
        File? input_bam
        File? input_bam_index
        String gvcf_basename
        String sampleName
        String? sampleGroup
        String? scatterName
        String? scatterIntervals
        File? FinalCallset
        File? FinalCallset_index
        # runtime
        String container
        Int? runtime_set_preemptible_tries
        Int? runtime_set_cpu
        Int? runtime_set_memory
        Int? runtime_set_disk
        Int? runtime_set_max_retries
        Boolean use_ssd = false
    }
    Float size_input_files = size(ref, "GB") + size(ref_dict, "GB") + size(input_bam, "GB") + size(input_bam_index, "GB") + size(FinalCallset, "GB") + size(FinalCallset_index, "GB")
    Int runtime_calculated_disk = ceil(size_input_files * 1.5) + 20
    Int command_mem_gb = select_first([runtime_set_memory, 8]) - 1
    String padding = if defined(FinalCallset) then "--interval-padding 100 \\" else " "
    command {
        gatk \
        HaplotypeCaller --java-options "-Xmx~{command_mem_gb}G" \
        -I ~{input_bam} \
        -R ~{ref} \
        -O ~{gvcf_basename}.g.vcf.gz \
        -L ~{scatterIntervals} ~{FinalCallset} \
        -ERC GVCF \
        -G StandardAnnotation \
        -G AS_StandardAnnotation \
        -G StandardHCAnnotation \
        -GQB 10 -GQB 20 -GQB 30 -GQB 40 -GQB 50 -GQB 60 -GQB 70 -GQB 80 -GQB 90 \
        ~{padding}
    }
    output {
        File output_gvcf = "~{gvcf_basename}" + ".g.vcf.gz"
        File output_gvcf_index = "~{gvcf_basename}" + ".g.vcf.gz.tbi"
        String output_gvcf_name = "~{sampleName}"
        String output_gvcf_scattername = "~{scatterName}"
        String output_gvcf_group = "~{sampleGroup}"
    }
    runtime {
        container: container
        cpu: select_first([runtime_set_cpu, 1])
        gpu: false
        memory: select_first([runtime_set_memory, 8]) + " GB"
        disks: "local-disk " + select_first([runtime_set_disk, runtime_calculated_disk]) + if use_ssd then " SSD" else " HDD"
        maxRetries: select_first([runtime_set_max_retries, 0])
        preemptible: select_first([runtime_set_preemptible_tries, 5])
        returnCodes: 0
     }
}

task GenomicsDBImport {
    input {
        Array[File] input_gvcfs
        String database_name
        String scatterIntervals
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
    Float size_input_files = size(input_gvcfs, "GB")
    Int runtime_calculated_disk = ceil(size_input_files * 2.5) + 20
    Int command_mem_gb = select_first([runtime_set_memory, 8]) - 1
    #Int merge_contigs_value = if defined(merge_contigs_into_num_partitions) then merge_contigs_into_num_partitions else "0"
    command <<<
    set -euo pipefail

        gatk IndexFeatureFile -I ~{sep=' && gatk IndexFeatureFile -I ' input_gvcfs}

        gatk \
        GenomicsDBImport --java-options "-Xmx~{command_mem_gb}G" \
        -V ~{sep=' -V ' input_gvcfs} \
        --genomicsdb-workspace-path ~{database_name} \
        -L ~{scatterIntervals}

        tar -cf ~{database_name}.tar ~{database_name}
    >>>
    output {
        File output_genomicsdb = "~{database_name}.tar"
    }
    runtime {
        container: container
        cpu: select_first([runtime_set_cpu, 1])
        gpu: false
        memory: select_first([runtime_set_memory, 8]) + " GB"
        disks: "local-disk " + select_first([runtime_set_disk, runtime_calculated_disk]) + if use_ssd then " SSD" else " HDD"
        maxRetries: select_first([runtime_set_max_retries, 0])
        preemptible: select_first([runtime_set_preemptible_tries, 5])
        returnCodes: 0
     }
}

task GenomicsDBImportFinal {
    input {
        Array[File] input_gvcfs
        Array[File] input_gvcfs_indexes
        String database_name
        File FinalCallset
        File FinalCallset_index
        # runtime
        String container
        Int? runtime_set_preemptible_tries
        Int? runtime_set_cpu
        Int? runtime_set_memory
        Int? runtime_set_disk
        Int? runtime_set_max_retries
        Boolean use_ssd = false
        Boolean use_ssd = false
    }
    Float size_input_files = size(input_gvcfs, "GB")
    Int runtime_calculated_disk = ceil(size_input_files * 2.5) + 20
    Int command_mem_gb = select_first([runtime_set_memory, 8]) - 1
    command <<<
    set -euo pipefail

        gatk \
        GenomicsDBImport --java-options "-Xmx~{command_mem_gb}G" \
        -V ~{sep=' -V ' input_gvcfs} \
        --genomicsdb-workspace-path ~{database_name} \
        -L ~{FinalCallset} \
        --interval-padding 100

        tar -cf ~{database_name}.tar ~{database_name}
    >>>
    output {
        File output_genomicsdb = "~{database_name}.tar"
    }
    runtime {
        container: container
        cpu: select_first([runtime_set_cpu, 1])
        gpu: false
        memory: select_first([runtime_set_memory, 8]) + " GB"
        disks: "local-disk " + select_first([runtime_set_disk, runtime_calculated_disk]) + if use_ssd then " SSD" else " HDD"
        maxRetries: select_first([runtime_set_max_retries, 0])
        preemptible: select_first([runtime_set_preemptible_tries, 5])
        returnCodes: 0
     }
}


task GenotypeGenomicsDB {
    input {
        File ref
        File ref_dict
        Array[File]+ ref_idxs
        String database_name
        File input_genomicsdb
        File? FinalCallset
        File? FinalCallset_index
        String? scatterIntervals
        String? groupName
        # runtime
        String container
        Int? runtime_set_preemptible_tries
        Int? runtime_set_cpu
        Int? runtime_set_memory
        Int? runtime_set_disk
        Int? runtime_set_max_retries
        Boolean use_ssd = false
    }
    Float size_input_files = size(ref, "GB") + size(ref_dict, "GB") + size(ref_idxs, "GB") + size(input_genomicsdb, "GB") + size(FinalCallset, "GB") + size(FinalCallset_index, "GB")
    Int runtime_calculated_disk = ceil(size_input_files * 2.5) + 20
    Int command_mem_gb = select_first([runtime_set_memory, 8]) - 1
    String final_options = if defined(FinalCallset) then "--include-non-variant-sites\\" else " "
    command <<<
    set -euo pipefail

        tar -xf ~{input_genomicsdb}

        gatk \
        GenotypeGVCFs --java-options "-Xmx~{command_mem_gb}G" \
        -R ~{ref} \
        -V gendb://~{database_name} \
        -L ~{scatterIntervals} ~{FinalCallset} \
        -O ~{database_name}.vcf.gz \
        -G StandardAnnotation -G AS_StandardAnnotation \
        ~{final_options}
        
    >>>
    output {
        String output_groupName = "~{groupName}"
        File output_genotypes = "~{database_name}.vcf.gz"
        File output_genotypes_index = "~{database_name}.vcf.gz.tbi"
    }
    runtime {
        container: container
        cpu: select_first([runtime_set_cpu, 1])
        gpu: false
        memory: select_first([runtime_set_memory, 8]) + " GB"
        disks: "local-disk " + select_first([runtime_set_disk, runtime_calculated_disk]) + if use_ssd then " SSD" else " HDD"
        maxRetries: select_first([runtime_set_max_retries, 0])
        preemptible: select_first([runtime_set_preemptible_tries, 5])
        returnCodes: 0
     }
}
