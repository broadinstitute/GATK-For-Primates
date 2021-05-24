version development

task ProduceFinalCallset {
    input {
        Array[File]? SNPs_hard_filtered
        Array[File]? SNPs_hard_filtered_index
        Array[File]? INDELs_hard_filtered
        Array[File]? INDELs_hard_filtered_index
        File? SNPs_recalibrated
        File? SNPs_recalibrated_index
        File? INDELs_recalibrated
        File? INDELs_recalibrated_index
        File ref_dict
        # runtime
        String container
        Int? runtime_set_preemptible_tries
        Int? runtime_set_cpu
        Int? runtime_set_memory
        Int? runtime_set_disk
        Int? runtime_set_max_retries
        Boolean use_ssd = false  
    }
    String a = if defined(SNPs_hard_filtered) then " -I " else ""
    String b = if defined(SNPs_recalibrated) then " -I " else ""
    String c = if defined(INDELs_hard_filtered) then " -I " else ""
    String d = if defined(INDELs_recalibrated) then " -I " else ""
    Float size_input_files = size(SNPs_hard_filtered, "GB") + size(SNPs_hard_filtered_index, "GB") + size(INDELs_hard_filtered, "GB") + size(INDELs_hard_filtered_index, "GB") + size(SNPs_recalibrated, "GB") + size(SNPs_recalibrated_index, "GB") + size(INDELs_recalibrated, "GB") + size(INDELs_recalibrated_index, "GB") + size(ref_dict, "GB")
    Int runtime_calculated_disk = ceil(size_input_files * 2) + 20
    Int command_mem_gb = runtime_set_memory - 1
    command <<<
        set -euo pipefail

        gatk MergeVcfs --java-options "-Xmx~{command_mem_gb}G" \
        ~{a} ~{sep=" -I " SNPs_hard_filtered} ~{b} ~{sep=" -I " SNPs_recalibrated} ~{c} ~{sep=" -I " INDELs_hard_filtered} ~{d} ~{sep=" -I " INDELs_recalibrated} \
        -O FinalCallset.vcf.gz

    >>>
    output {
        File FinalCallset = "FinalCallset.vcf.gz"
        File FinalCallset_index = "FinalCallset.vcf.gz.tbi"
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
