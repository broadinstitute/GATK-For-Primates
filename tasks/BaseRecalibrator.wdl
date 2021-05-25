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
        # runtime
        String container
        Int? runtime_set_preemptible_tries
        Int? runtime_set_cpu
        Int? runtime_set_memory
        Int? runtime_set_disk
        Int? runtime_set_max_retries
        Boolean use_ssd = false
    }
    Float size_input_files = size(input_bam, "GB") + size(input_bam_index, "GB") + size(ref, "GB") + size(ref_dict, "GB") + size(ref_idxs, "GB")
    Int runtime_calculated_disk = ceil(size_input_files * 2.5) + 20
    Int command_mem_gb = select_first([runtime_set_memory, 8]) - 1
    command <<<
        set -euo pipefail
        
        #GENERATE BEFORE TABLE FROM THE UNRECALIBRATED BAM FILE USING THE HARD FILTERED CALLS FROM THE CURRENT/LATEST PASS
        gatk \
        BaseRecalibrator --java-options "-Xmx~{command_mem_gb}G" \
        -I ~{input_bam} \
        -R ~{ref} \
        --known-sites ~{sep=" --known-sites " input_SNP_sites} \
        --known-sites ~{sep=" --known-sites " input_INDEL_sites} \
        -O ~{sampleName}.before.table

        #APPLY BQSR TO ADJUST THE QUALITY SCORES
        gatk \
        ApplyBQSR --java-options "-Xmx~{command_mem_gb}G" \
        -I ~{input_bam} \
        -R ~{ref} \
        --bqsr-recal-file ~{sampleName}.before.table \
        -O ~{sampleName}_recalibrated.bam

        #GENERATE AFTER TABLE
        gatk \
        BaseRecalibrator --java-options "-Xmx~{command_mem_gb}G" \
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

task AnalyzeCovariates {
    input {
        String sampleName
        File table_before
        File table_after
        # runtime
        String container
        Int? runtime_set_preemptible_tries
        Int? runtime_set_cpu
        Int? runtime_set_memory
        Int? runtime_set_disk
        Int? runtime_set_max_retries
        Boolean use_ssd = false
    }
    Int command_mem_gb = select_first([runtime_set_memory, 8]) - 1
    command <<<
        set -euo pipefail

        # The required R packages are not in Broad's production GITC image
        R --vanilla << CODE
        install.packages("gplots", repos="http://cran.us.r-project.org")
        install.packages("gsalib", repos="http://cran.us.r-project.org")
        install.packages("reshape", repos="http://cran.us.r-project.org")
        CODE
        
        gatk --java-options "-Xmx~{command_mem_gb}G" \
        AnalyzeCovariates \
        -before ~{table_before} \
        -after ~{table_after} \
        -plots ~{sampleName}.AnalyzeCovariates_plots.pdf \
        --use-jdk-deflater \
        --use-jdk-inflater
    >>>
    output {
        File plots = "~{sampleName}.AnalyzeCovariates_plots.pdf"
    }
    runtime {
        container: container
        cpu: select_first([runtime_set_cpu, 1])
        gpu: false
        memory: select_first([runtime_set_memory, 8]) + " GB"
        disks: "local-disk " + select_first([runtime_set_disk, 25]) + if use_ssd then " SSD" else " HDD"
        maxRetries: select_first([runtime_set_max_retries, 0])
        preemptible: select_first([runtime_set_preemptible_tries, 5])
        returnCodes: 0
     }
}
