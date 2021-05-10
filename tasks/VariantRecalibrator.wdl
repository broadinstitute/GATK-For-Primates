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
        String docker_image
    }
    command <<<
        set -euo pipefail

        if (type == "SNP") {

            gatk \
            VariantRecalibrator \
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
            ApplyVQSR \
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
            VariantRecalibrator \
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
            ApplyVQSR \
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

        gatk MakeSitesOnlyVcf \
        -I Cohort_~{type}.recalibrated.vcf.gz \
        -O Cohort_~{type}.recalibrated_sites_only.vcf.gz


    >>>
    output {
        File output_recalibrated_sites_only = "Cohort_~{type}.recalibrated_sites_only.vcf.gz"
        File output_recalibrated_sites_only_index = "Cohort_~{type}.recalibrated_sites_only.vcf.gz.tbi"
    }
}