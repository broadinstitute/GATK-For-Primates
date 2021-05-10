version development

task BaseRecalibrator {
    input {
        File ref
        File ref_dict
        Array[File]+ ref_idxs
        File input_bam
        File input_bam_index
        String sampleName
        Array[File] input_SNP_sites
        Array[File] input_SNP_sites_indexes
        Array[File] input_INDEL_sites
        Array[File] input_INDEL_sites_indexes
        String docker_image
    }
    command <<<
        set -euo pipefail
        
        #GENERATE BEFORE TABLE FROM THE UNRECALIBRATED BAM FILE USING THE HARD FILTERED CALLS FROM THE CURRENT/LATEST PASS
        gatk \
        BaseRecalibrator \
        -I ~{input_bam} \
        -R ~{ref} \
        --known-sites ~{sep=" --known-sites " input_SNP_sites} \
        --known-sites ~{sep=" --known-sites " input_INDEL_sites} \
        -O ~{sampleName}.before.table

        #APPLY BQSR TO ADJUST THE QUALITY SCORES
        gatk \
        ApplyBQSR \
        -I ~{input_bam} \
        -R ~{ref} \
        --bqsr-recal-file ~{sampleName}.before.table \
        -O ~{sampleName}_recalibrated.bam

        #GENERATE AFTER TABLE
        gatk \
        BaseRecalibrator \
        -I ~{sampleName}_recalibrated.bam \
        -R ~{ref} \
        --known-sites ~{sep=" --known-sites " input_SNP_sites} \
        --known-sites ~{sep=" --known-sites " input_INDEL_sites} \
        -O ~{sampleName}.after.table

        #GENERATE PLOTS
        #gatk \
        #AnalyzeCovariates \
        #-before ~{sampleName}.before.table \
        #-after ~{sampleName}.after.table \
        #-plots ~{sampleName}.AnalyzeCovariates_plots.pdf \
        #--use-jdk-deflater \
        #--use-jdk-inflater
    >>>
    output {
        File table_before = "~{sampleName}.before.table"
        File table_after = "~{sampleName}.after.table"
        #File plots = "~{sampleName}.AnalyzeCovariates_plots.pdf"
        File recalibrated_bam = "~{sampleName}_recalibrated.bam"
        File recalibrated_bam_index = "~{sampleName}_recalibrated.bai"
    }
}
