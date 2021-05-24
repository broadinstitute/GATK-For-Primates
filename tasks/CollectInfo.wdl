version development

task collectScatters {
    input {
        String scatterName
        String scatterIntervals
        # runtime
        String container
        Int? runtime_set_preemptible_tries
        Int? runtime_set_cpu
        Int? runtime_set_memory
        Int? runtime_set_disk
        Int? runtime_set_max_retries
        Boolean use_ssd = false
    }
    command {
    }
    output {
        String allScatterNames = scatterName
        String allScatterIntervals = scatterIntervals
    }
    runtime {
        container: container
        cpu: select_first([runtime_set_cpu, 1])
        gpu: false
        memory: select_first([runtime_set_memory, 5]) + " GB"
        disks: "local-disk " + select_first([runtime_set_disk, 10]) + if use_ssd then " SSD" else " HDD"
        maxRetries: select_first([runtime_set_max_retries, 0])
        preemptible: select_first([runtime_set_preemptible_tries, 5])
        returnCodes: 0
     }
}

task collectAllGroupNames {
    input {
        String groupName
        String sampleName
        # runtime
        String container
        Int? runtime_set_preemptible_tries
        Int? runtime_set_cpu
        Int? runtime_set_memory
        Int? runtime_set_disk
        Int? runtime_set_max_retries
        Boolean use_ssd = false
    }
    command {
    }
    output {
        String allGroupNames = groupName
        String allSampleNames = sampleName
    }
    runtime {
        container: container
        cpu: select_first([runtime_set_cpu, 1])
        gpu: false
        memory: select_first([runtime_set_memory, 5]) + " GB"
        disks: "local-disk " + select_first([runtime_set_disk, 10]) + if use_ssd then " SSD" else " HDD"
        maxRetries: select_first([runtime_set_max_retries, 0])
        preemptible: select_first([runtime_set_preemptible_tries, 5])
        returnCodes: 0
     }
}


task collectExistingBamInfo {
    input {
        String sampleName
        String sampleGroup
        File? sampleBAM
        File? sampleBAI
        # runtime
        String container
        Int? runtime_set_preemptible_tries
        Int? runtime_set_cpu
        Int? runtime_set_memory
        Int? runtime_set_disk
        Int? runtime_set_max_retries
        Boolean use_ssd = false
    }
    command <<<
        ## The only purpose of this is to coerce the ? inputs to outputs
    >>>
    output {
        String allNames = sampleName
        String allGroups = sampleGroup
        String allBams = sampleBAM
        String allBais = sampleBAI
    }
    runtime {
        container: container
        cpu: select_first([runtime_set_cpu, 1])
        gpu: false
        memory: select_first([runtime_set_memory, 5]) + " GB"
        disks: "local-disk " + select_first([runtime_set_disk, 10]) + if use_ssd then " SSD" else " HDD"
        maxRetries: select_first([runtime_set_max_retries, 0])
        preemptible: select_first([runtime_set_preemptible_tries, 5])
        returnCodes: 0
     }
}

task coerceFile {
    input {
        File? input_file
        # runtime
        String container
        Int? runtime_set_preemptible_tries
        Int? runtime_set_cpu
        Int? runtime_set_memory
        Int? runtime_set_disk
        Int? runtime_set_max_retries
        Boolean use_ssd = false
    }
    command <<<
        ## The only purpose of this task is to coerce File? into File
    >>>
    output {
        File coerced = input_file
    }
    runtime {
        container: container
        cpu: select_first([runtime_set_cpu, 1])
        gpu: false
        memory: select_first([runtime_set_memory, 5]) + " GB"
        disks: "local-disk " + select_first([runtime_set_disk, 10]) + if use_ssd then " SSD" else " HDD"
        maxRetries: select_first([runtime_set_max_retries, 0])
        preemptible: select_first([runtime_set_preemptible_tries, 5])
        returnCodes: 0
     }
}
