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
        Int? machine_mem_gb
        Int? disk_space_gb
        Int? preemptible_tries
        Boolean use_ssd = false
    }
    Float size_input_files = size(input_bam, "GB") + size(input_bam_index, "GB") + size(ref, "GB") + size(ref_dict, "GB") + size(ref_idxs, "GB")
    Int disk_size = ceil(size_input_files * 2.5) + 20
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
    runtime {
        docker: docker_image
        memory: select_first([machine_mem_gb, 8]) + " GB"
        disks: "local-disk " + select_first([disk_space_gb, disk_size]) + if use_ssd then " SSD" else " HDD"
        preemptible: select_first([preemptible_tries, 5])
    }
}

task AnalyzeCovariates {
    input {
        String sampleName
        File table_before
        File table_after
        String docker_image
        Int? machine_mem_gb
        Int? disk_space_gb
        Int? preemptible_tries
        Boolean use_ssd = false
    }
    command <<<
        set -euo pipefail

        gatk \
        AnalyzeCovariates \
        -before ~{sampleName}.before.table \
        -after ~{sampleName}.after.table \
        -plots ~{sampleName}.AnalyzeCovariates_plots.pdf \
        --use-jdk-deflater \
        --use-jdk-inflater
    >>>
    output {
        File plots = "~{sampleName}.AnalyzeCovariates_plots.pdf"
    }
    runtime {
        docker: docker_image
        memory: select_first([machine_mem_gb, 8]) + " GB"
        disks: "local-disk " + select_first([disk_space_gb, 25]) + if use_ssd then " SSD" else " HDD"
        preemptible: select_first([preemptible_tries, 5])
    }
}
