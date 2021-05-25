version development

task FastqToBwaMem {
    input {
        File ref
        File ref_dict
        Array[File]+ ref_idxs
        String sampleName
        String? rgID
        String? rgSM
        String? rgLB
        File? R1
        File? R2
        String execute_bwa
        # runtime
        String container
        Int? runtime_set_preemptible_tries
        Int? runtime_set_cpu
        Int? runtime_set_memory
        Int? runtime_set_disk
        Int? runtime_set_max_retries
        Boolean use_ssd = false
    }
    Float size_input_files = size(R1, "GB") + size(R2, "GB") + size(ref, "GB") + size(ref_dict, "GB") + size(ref_idxs, "GB")
    Int runtime_calculated_disk = ceil(size_input_files * 1.8) + 20
    command {
        ${execute_bwa} mem -K 100000000 -Y -v 3 -R "@RG\tID:${rgID}\tPL:ILLUMINA\tLB:${rgLB}\tSM:${rgSM}" ${ref} ${R1} ${R2} | samtools sort -n -o ~{sampleName}.bam -
    }
    output {
        File output_bam = "~{sampleName}.bam"
    }
    runtime {
        container: container
        cpu: select_first([runtime_set_cpu, 1])
        gpu: false
        memory: select_first([runtime_set_memory, 16]) + " GB"
        disks: "local-disk " + select_first([runtime_set_disk, runtime_calculated_disk]) + if use_ssd then " SSD" else " HDD"
        maxRetries: select_first([runtime_set_max_retries, 0])
        preemptible: select_first([runtime_set_preemptible_tries, 5])
        returnCodes: 0
     }
}

task MarkDuplicatesSpark {
    input {
        File input_bam
        String sampleName
        Boolean flowcell_patterned
        # runtime
        String container
        Int? runtime_set_preemptible_tries
        Int? runtime_set_cpu
        Int? runtime_set_memory
        Int? runtime_set_disk
        Int? runtime_set_max_retries
        Boolean use_ssd = false
    }
    Float size_input_files = size(input_bam, "GB")
    Int runtime_calculated_disk = ceil(size_input_files * 2.5) + 20
    String pixel_distance = if flowcell_patterned then "2500" else "100"
    Int command_mem_gb = select_first([runtime_set_memory, 24]) - 6
    command {
        gatk \
        MarkDuplicatesSpark --java-options "-Xmx~{command_mem_gb}G" \
        -I ~{input_bam} \
        --optical-duplicate-pixel-distance ~{pixel_distance} \
        -O ~{sampleName}.dedup.bam \
        -M ~{sampleName}.metrics.txt
    }
    output {
        File output_bam = "~{sampleName}.dedup.bam"
        File duplication_metrics = "~{sampleName}.metrics.txt"
    }
    runtime {
        container: container
        cpu: select_first([runtime_set_cpu, 1])
        gpu: false
        memory: select_first([runtime_set_memory, 24]) + " GB"
        disks: "local-disk " + select_first([runtime_set_disk, runtime_calculated_disk]) + if use_ssd then " SSD" else " HDD"
        maxRetries: select_first([runtime_set_max_retries, 0])
        preemptible: select_first([runtime_set_preemptible_tries, 5])
        returnCodes: 0
     }
}

task SortAndFixTags {
    input {
        File ref
        File ref_dict
        Array[File]+ ref_idxs
        File input_bam
        String sampleName
        String sampleGroup
        # runtime
        String container
        Int? runtime_set_preemptible_tries
        Int? runtime_set_cpu
        Int? runtime_set_memory
        Int? runtime_set_disk
        Int? runtime_set_max_retries
        Boolean use_ssd = false
    }
    Float size_input_files = size(input_bam, "GB") + size(ref, "GB") + size(ref_dict, "GB") + size(ref_idxs, "GB")
    Int runtime_calculated_disk = ceil(size_input_files * 2) + 20
    Int command_mem_gb = select_first([runtime_set_memory, 8]) - 1
    command {
        gatk \
        SetNmMdAndUqTags --java-options "-Xmx~{command_mem_gb}G" \
        -I ~{input_bam} \
        -R ~{ref} \
        -O ~{sampleName}.dedup.tagged.bam \
        --CREATE_INDEX true
    }
    output {
        File output_bam = "~{sampleName}.dedup.tagged.bam"
        File output_bam_index = "~{sampleName}.dedup.tagged.bai"
        String output_name = "~{sampleName}"        
        String output_group = "~{sampleGroup}"
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
