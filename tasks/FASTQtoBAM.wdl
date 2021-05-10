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
        String docker_image
        Int? machine_mem_gb
        Int? disk_space_gb
        Int? preemptible_tries
        Boolean use_ssd = false
    }
    Float size_input_files = size(R1, "GB") + size(R2, "GB") + size(ref, "GB") + size(ref_dict, "GB") + size(ref_idxs, "GB")
    Int disk_size = ceil(size_input_files * 1.8) + 20
    command {
        ${execute_bwa} mem -K 100000000 -Y -v 3 -R "@RG\tID:${rgID}\tPL:ILLUMINA\tLB:${rgLB}\tSM:${rgSM}" ${ref} ${R1} ${R2} | samtools sort -n -o ~{sampleName}.bam -
    }
    output {
        File output_bam = "~{sampleName}.bam"
    }
    runtime {
        docker: docker_image
        memory: select_first([machine_mem_gb, 16]) + " GB"
        disks: "local-disk " + select_first([disk_space_gb, disk_size]) + if use_ssd then " SSD" else " HDD"
        preemptible: select_first([preemptible_tries, 5])
     }
}

task MarkDuplicatesSpark {
    input {
        File input_bam
        String sampleName
        Boolean flowcell_patterned
        String docker_image
        Int? machine_mem_gb
        Int? disk_space_gb
        Int? preemptible_tries
        Boolean use_ssd = false
    }
    Float size_input_files = size(input_bam, "GB")
    Int disk_size = ceil(size_input_files * 2.5) + 20
    String pixel_distance = if flowcell_patterned then "2500" else "100"
        command {
        gatk \
        MarkDuplicatesSpark \
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
        docker: docker_image
        memory: select_first([machine_mem_gb, 16]) + " GB"
        disks: "local-disk " + select_first([disk_space_gb, disk_size]) + if use_ssd then " SSD" else " HDD"
        preemptible: select_first([preemptible_tries, 5])
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
        String docker_image
        Int? machine_mem_gb
        Int? disk_space_gb
        Int? preemptible_tries
        Boolean use_ssd = false
    }
    Float size_input_files = size(input_bam, "GB") + size(ref, "GB") + size(ref_dict, "GB") + size(ref_idxs, "GB")
    Int disk_size = ceil(size_input_files * 2) + 20
    command {
        gatk \
        SetNmMdAndUqTags \
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
        docker: docker_image
        memory: select_first([machine_mem_gb, 8]) + " GB"
        disks: "local-disk " + select_first([disk_space_gb, disk_size]) + if use_ssd then " SSD" else " HDD"
        preemptible: select_first([preemptible_tries, 5])
    }
}
