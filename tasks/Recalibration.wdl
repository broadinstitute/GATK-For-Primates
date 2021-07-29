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
## *** TASK: baseRecalibrator ***
##########################################################################
## Generates before table, applies BQSR, then generates after table.
##########################################################################

task baseRecalibrator {
    input {
        File ref
        File ref_dict
        File ref_fai
        File input_bam
        File input_bam_index
        Boolean cram_not_bam
        String sampleName
        File high_confidence_sites_vcf
        File high_confidence_sites_vcf_index
        # runtime
        String container
        Int? runtime_set_preemptible_tries
        Int? runtime_set_cpu
        Int? runtime_set_memory
        Int? runtime_set_disk
        Int? runtime_set_max_retries
        Boolean use_ssd = false
    }
    Float size_input_files = size(input_bam, "GB") + size(input_bam_index, "GB") + size(ref, "GB") + size(ref_dict, "GB") + size(ref_fai, "GB") + size(high_confidence_sites_vcf, "GB") + size(high_confidence_sites_vcf_index, "GB")
    Int runtime_calculated_disk = ceil(size_input_files * 2.5)
    Int command_mem_gb = select_first([runtime_set_memory, 6]) - 2
    String extension = if cram_not_bam then "cram" else "bam"
    String extension_index = if cram_not_bam then "crai" else "bai"
    command <<<
        set -euo pipefail
        
        #GENERATE BEFORE TABLE FROM THE UNRECALIBRATED BAM FILE USING THE HARD FILTERED CALLS FROM THE CURRENT/LATEST PASS
        gatk \
        BaseRecalibrator --java-options "-Xmx~{command_mem_gb}G" \
        -I ~{input_bam} \
        -R ~{ref} \
        --known-sites ~{high_confidence_sites_vcf} \
        -O ~{sampleName}.before.table

        #APPLY BQSR TO ADJUST THE QUALITY SCORES
        gatk \
        ApplyBQSR --java-options "-Xmx~{command_mem_gb}G" \
        -I ~{input_bam} \
        -R ~{ref} \
        --bqsr-recal-file ~{sampleName}.before.table \
        -O ~{sampleName}_recalibrated.~{extension}

        samtools index ~{sampleName}_recalibrated.~{extension} ~{sampleName}_recalibrated.~{extension_index}

        #GENERATE AFTER TABLE
        gatk \
        BaseRecalibrator --java-options "-Xmx~{command_mem_gb}G" \
        -I ~{sampleName}_recalibrated.~{extension} \
        -R ~{ref} \
        --known-sites ~{high_confidence_sites_vcf} \
        -O ~{sampleName}.after.table

    >>>
    output {
        File table_before = "~{sampleName}.before.table"
        File table_after = "~{sampleName}.after.table"
        File recalibrated_bam = "~{sampleName}_recalibrated.~{extension}"
        File recalibrated_bam_index = "~{sampleName}_recalibrated.~{extension_index}"
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
## *** TASK: analyzeCovariates ***
##########################################################################
## Plots the covariates.
##########################################################################

task analyzeCovariates {
    input {
        String sampleName
        File table_before
        File table_after
        # runtime
        String container
        String path_to_gitc_gatk
        Int? runtime_set_preemptible_tries
        Int? runtime_set_cpu
        Int? runtime_set_memory
        Int? runtime_set_disk
        Int? runtime_set_max_retries
        Boolean use_ssd = false
    }
    Int command_mem_gb = select_first([runtime_set_memory, 6]) - 1
    command <<<
        set -euo pipefail

        # The required R packages are not in Broad's production GITC image
        R --vanilla << CODE
        install.packages("gplots", repos="http://cran.us.r-project.org")
        install.packages("gsalib", repos="http://cran.us.r-project.org")
        install.packages("reshape", repos="http://cran.us.r-project.org")
        CODE
        
        ~{path_to_gitc_gatk}gatk --java-options "-Xmx~{command_mem_gb}G" \
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
        docker: container
        cpu: select_first([runtime_set_cpu, 1])
        gpu: false
        memory: select_first([runtime_set_memory, 6]) + " GB"
        disks: "local-disk " + select_first([runtime_set_disk, 3]) + if use_ssd then " SSD" else " HDD"
        maxRetries: select_first([runtime_set_max_retries, 0])
        preemptible: select_first([runtime_set_preemptible_tries, 5])
        returnCodes: 0
     }
}


##########################################################################
## *** TASK: variantRecalibrator ***
##########################################################################
## Recalibrates variants according to VQSR.
##########################################################################

task variantRecalibrator {
    input {
        File ref
        File ref_dict
        File ref_fai
        String groupName
        String type
        File input_vcf_sites_only
        File input_vcf_sites_only_index
        File? truth_set
        File? truth_set_index
        File training_set
        File training_set_index
        # Runtime
        String container
        Int? runtime_set_preemptible_tries
        Int? runtime_set_cpu
        Int? runtime_set_memory
        Int? runtime_set_disk
        Int? runtime_set_max_retries
        Boolean use_ssd = false
    }
    Float size_input_files = size(ref, "GB") + size(ref_dict, "GB") + size(ref_fai, "GB") + size(truth_set, "GB") + size(truth_set_index, "GB") + size(training_set, "GB") + size(training_set_index, "GB") + size(input_vcf_sites_only, "GB") + size(input_vcf_sites_only_index, "GB")
    Int runtime_calculated_disk = ceil(size_input_files * 2) + 20
    Int command_mem_gb = select_first([runtime_set_memory, 8]) - 2
    Int max_gaussians = if "~{type}" == "INDEL" then 4 else 6
    command <<<
        set -euo pipefail

            gatk \
            VariantRecalibrator --java-options "-Xmx~{command_mem_gb}G" \
            -R ~{ref} \
            -V ~{input_vcf_sites_only} \
            -O finalCallset_~{groupName}_~{type}s.recal \
            --trust-all-polymorphic \
            -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 97.0 -tranche 96.0 -tranche 95.0 -tranche 94.0 -tranche 93.5 -tranche 93.0 -tranche 92.0 -tranche 91.0 -tranche 90.0 \
            -AS \
            -mode ~{type} \
            --max-gaussians ~{max_gaussians} \            
            -resource:truth_set_~{type}s,known=false,training=true,truth=true,prior=12.0 ~{truth_set} \
            -resource:hard_filtered_~{type}s,known=false,training=true,truth=false,prior=10.0 ~{training_set} \
            --tranches-file finalCallset_~{groupName}_~{type}s.AS.tranches \
            --rscript-file finalCallset_~{groupName}_~{type}s.plots.AS.R

    >>>
    output {
        File output_recalibrated_recal = "finalCallset_~{groupName}_~{type}s.recal"
        File output_recalibrated_tranches = "finalCallset_~{groupName}_~{type}s.AS.tranches"
        File output_recalibrated_rscript = "finalCallset_~{groupName}_~{type}s.plots.AS.R"
    }
    runtime {
        docker: container
        cpu: select_first([runtime_set_cpu, 1])
        gpu: false
        memory: select_first([runtime_set_memory, 8]) + " GB"
        disks: "local-disk " + select_first([runtime_set_disk, runtime_calculated_disk]) + if use_ssd then " SSD" else " HDD"
        maxRetries: select_first([runtime_set_max_retries, 0])
        preemptible: select_first([runtime_set_preemptible_tries, 1])
        returnCodes: 0
     }
}

##########################################################################
## *** TASK: applyVQSR ***
##########################################################################
## Applies VQSR.
##########################################################################

task applyVQSR {
    input {
        File ref
        File ref_dict
        File ref_fai
        String groupName
        String type
        File input_vcf
        File input_vcf_index
        File recalibrated_recal
        File recalibrated_tranches
        File recalibrated_rscript
        Boolean exclude_filtered
        # Runtime
        String container
        Int? runtime_set_preemptible_tries
        Int? runtime_set_cpu
        Int? runtime_set_memory
        Int? runtime_set_disk
        Int? runtime_set_max_retries
        Boolean use_ssd = false
    }
    Float size_input_files = size(ref, "GB") + size(ref_dict, "GB") + size(input_vcf, "GB") + size(input_vcf_index, "GB") + size(recalibrated_recal, "GB") + size(recalibrated_tranches, "GB") + size(recalibrated_rscript, "GB")
    Int runtime_calculated_disk = ceil(size_input_files * 2.5)
    Int command_mem_gb = select_first([runtime_set_memory, 7]) - 2
    String exclude_filtered_param = if exclude_filtered then "--exclude-filtered" else ""
    command <<<
        set -euo pipefail

            gatk \
            ApplyVQSR --java-options "-Xmx~{command_mem_gb}G" \
            -R ~{ref} \
            -V ~{input_vcf} \
            -O finalCallset_~{groupName}_~{type}s.recalibrated.vcf.gz
            -AS \
            --recal-file finalCallset.~{type}s.recal \
            --truth-sensitivity-filter-level 99.0 \
            --tranches-file finalCallset.~{type}s.AS.tranches \
            --create-output-variant-index true \
            -mode ~{type} \
            ~{exclude_filtered}

    >>>
    output {
        File output_recalibrated_vcf = "Cohort_~{groupName}.recalibrated.vcf.gz"
        File output_recalibrated_vcf_index = "Cohort_~{groupName}.recalibrated.vcf.gz.tbi"
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
