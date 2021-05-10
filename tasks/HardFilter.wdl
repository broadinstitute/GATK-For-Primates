version development

## Need to add AS annotation flags here
task SNPs {
    input {
        File ref
        File ref_dict
        Array[File]+ ref_idxs
        File input_genotypes
        File input_genotypes_index
        String groupName
        String scatterName
        String scatterIntervals
        String docker_image
    }
    command {

        gatk \
        SelectVariants \
        -R ~{ref} \
        -V ~{input_genotypes} \
        --select-type-to-include SNP \
        -O ~{groupName}_~{scatterName}_SNPs_unfiltered.vcf.gz

        gatk \
        VariantFiltration \
        -R ~{ref} \
        -V ~{groupName}_~{scatterName}_SNPs_unfiltered.vcf.gz \
        -O ~{groupName}_~{scatterName}_SNPs_unfiltered_and_filtered.vcf.gz \
        -filter "QD < 2.0" --filter-name "QD2" \
        -filter "QUAL < 30.0" --filter-name "QUAL30" \
        -filter "SOR > 3.0" --filter-name "SOR3" \
        -filter "FS > 60.0" --filter-name "FS60" \
        -filter "MQ < 40.0" --filter-name "MQ40" \
        -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
        -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
        --verbosity ERROR

        gatk \
        SelectVariants \
        -R ~{ref} \
        -V ~{groupName}_~{scatterName}_SNPs_unfiltered_and_filtered.vcf.gz \
        -O ~{groupName}_~{scatterName}_SNPs_filtered.vcf.gz \
        --exclude-filtered
        
    }
    output {
        File output_filtered_SNPs = "~{groupName}_~{scatterName}_SNPs_filtered.vcf.gz"
        File output_filtered_SNPs_index = "~{groupName}_~{scatterName}_SNPs_filtered.vcf.gz.tbi"
        String output_groupName = "~{groupName}"
    }
}

## Need to add AS annotation flags here
task INDELs {
    input {
        File ref
        File ref_dict
        Array[File]+ ref_idxs
        File input_genotypes
        File input_genotypes_index
        String groupName
        String scatterName
        String scatterIntervals
        String docker_image
    }
    command {

        gatk \
        SelectVariants \
        -R ~{ref} \
        -V ~{input_genotypes} \
        --select-type-to-include MIXED \
        -O ~{groupName}_~{scatterName}_INDELs_unfiltered.vcf.gz

        gatk \
        VariantFiltration \
        -R ~{ref} \
        -V ~{groupName}_~{scatterName}_INDELs_unfiltered.vcf.gz \
        -O ~{groupName}_~{scatterName}_INDELs_unfiltered_and_filtered.vcf.gz \
        -filter "QD < 2.0" --filter-name "QD2" \
        -filter "QUAL < 30.0" --filter-name "QUAL30" \
        -filter "FS > 200.0" --filter-name "FS200" \
        -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
        --verbosity ERROR

        gatk \
        SelectVariants \
        -R ~{ref} \
        -V ~{groupName}_~{scatterName}_INDELs_unfiltered_and_filtered.vcf.gz \
        -O ~{groupName}_~{scatterName}_INDELs_filtered.vcf.gz \
        --exclude-filtered
        
    }
    output {
        File output_filtered_INDELs = "~{groupName}_~{scatterName}_INDELs_filtered.vcf.gz"
        File output_filtered_INDELs_index = "~{groupName}_~{scatterName}_INDELs_filtered.vcf.gz.tbi"
        String output_groupName = "~{groupName}"
    }
}