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
    }
    command {
        ${execute_bwa} mem -K 100000000 -Y -v 3 -R "@RG\tID:${rgID}\tPL:ILLUMINA\tLB:${rgLB}\tSM:${rgSM}" ${ref} ${R1} ${R2} | samtools sort -n -o ~{sampleName}.bam -

        
    }
    output {
        File output_bam = "~{sampleName}.bam"
        
    }
    runtime {
        docker: docker_image
    }
}

task MarkDuplicatesSpark {
    input {
        String docker_image
        File input_bam
        String sampleName
        Boolean flowcell_patterned
    }
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
    #runtime {
    #    docker: docker_image
    #}
}

task SortAndFixTags {
    input {
        File ref
        File ref_dict
        Array[File]+ ref_idxs
        String docker_image
        File input_bam
        String sampleName
        String sampleGroup
    }
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
    }
}