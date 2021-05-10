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
    }
    String a = if defined(SNPs_hard_filtered) then " -I " else ""
    String b = if defined(SNPs_recalibrated) then " -I " else ""
    String c = if defined(INDELs_hard_filtered) then " -I " else ""
    String d = if defined(INDELs_recalibrated) then " -I " else ""
    command <<<
        set -euo pipefail

        gatk MergeVcfs \
        ~{a} ~{sep=" -I " SNPs_hard_filtered} ~{b} ~{sep=" -I " SNPs_recalibrated} ~{c} ~{sep=" -I " INDELs_hard_filtered} ~{d} ~{sep=" -I " INDELs_recalibrated} \
        -O FinalCallset.vcf.gz

    >>>
    output {
        File FinalCallset = "FinalCallset.vcf.gz"
        File FinalCallset_index = "FinalCallset.vcf.gz.tbi"
    }
}