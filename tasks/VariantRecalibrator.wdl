version development

task VariantRecalibrator {
    input {
        File ref
        File ref_dict
        Array[File]+ ref_idxs
        File? truth_set
        Array[File] input_sites
        Array[File] input_sites_indexes
        String type
        # runtime
        String container
        Int? runtime_set_preemptible_tries
        Int? runtime_set_cpu
        Int? runtime_set_memory
        Int? runtime_set_disk
        Int? runtime_set_max_retries
        Boolean use_ssd = false
    }
    Float size_input_files = size(ref, "GB") + size(ref_dict, "GB") + size(ref_idxs, "GB") + size(truth_set, "GB") + size(input_sites, "GB") + size(input_sites_indexes, "GB")
    Int runtime_calculated_disk = ceil(size_input_files * 2) + 20
    Int command_mem_gb = select_first([runtime_set_memory, 8]) - 1
    command <<<
        set -euo pipefail

        if (type == "SNP") {

            gatk \
            VariantRecalibrator --java-options "-Xmx~{command_mem_gb}G" \
            -V ~{sep=" -V " input_sites} \
            --trust-all-polymorphic \
            -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.8 -tranche 99.6 -tranche 99.5 -tranche 99.4 -tranche 99.3 -tranche 99.0 -tranche 98.0 -tranche 97.0 -tranche 90.0 \
            -an QD -an MQRankSum -an ReadPosRankSum -an FS -an MQ -an SOR -an DP \
            -mode SNP \
            --max-gaussians 6 \
            -resource:truth_set_snps,known=false,training=true,truth=true,prior=~{truth_set} \
            -resource:hard_filtered,known=false,training=true,truth=false,prior=~{input_sites} \
            -O cohort_snps.recal \
            --tranches-file cohort_snps.tranches \
            --use-allele-specific-annotations

            gatk \
            ApplyVQSR --java-options "-Xmx~{command_mem_gb}G" \
            -V ~{sep=" -V " input_sites} \
            --recal-file cohort_snps.recal \
            --tranches-file cohort_snps.tranches \
            --truth-sensitivity-filter-level 99.7 \
            --create-output-variant-index true \
            -mode SNP \
            --exclude-filtered \
            --use-allele-specific-annotations \
            -O Cohort_~{type}.recalibrated.vcf.gz

        }

        if (type == "INDEL") {

            gatk \
            VariantRecalibrator --java-options "-Xmx~{command_mem_gb}G" \
            -V ~{sep=" -V " input_sites} \
            --trust-all-polymorphic \
            -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 97.0 -tranche 96.0 -tranche 95.0 -tranche 94.0 -tranche 93.5 -tranche 93.0 -tranche 92.0 -tranche 91.0 -tranche 90.0 \
            -an FS -an ReadPosRankSum -an MQRankSum -an QD -an SOR -an DP \      
            -mode MIXED \
            --max-gaussians 4 \
            -resource:truth_set_indels,known=false,training=true,truth=true,prior=~{truth_set} \
            -resource:hard_filtered,known=false,training=true,truth=false,prior=~{input_sites} \
            -O cohort_indels.recal \
            --use-allele-specific-annotations \
            --tranches-file cohort_indels.tranches

            gatk \
            ApplyVQSR --java-options "-Xmx~{command_mem_gb}G" \
            -V ~{sep=" -V " input_sites} \
            --recal-file cohort_indels.recal \
            --tranches-file cohort_indels.tranches \
            --truth-sensitivity-filter-level 99.7 \
            --create-output-variant-index true \
            -mode MIXED \
            --exclude-filtered \
            --use-allele-specific-annotations \
            -O Cohort_~{type}.recalibrated.vcf.gz

        }

        gatk MakeSitesOnlyVcf --java-options "-Xmx~{command_mem_gb}G" \
        -I Cohort_~{type}.recalibrated.vcf.gz \
        -O Cohort_~{type}.recalibrated_sites_only.vcf.gz

    >>>
    output {
        File output_recalibrated_sites_only = "Cohort_~{type}.recalibrated_sites_only.vcf.gz"
        File output_recalibrated_sites_only_index = "Cohort_~{type}.recalibrated_sites_only.vcf.gz.tbi"
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
