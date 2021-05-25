version development

task FinalizeInputs {
    input {
        Array[Array[String]]? tsv_of_new_bams
        Array[Array[String]]? tsv_of_existing_bams
        # runtime
        String container
        Int? runtime_set_preemptible_tries
        Int? runtime_set_cpu
        Int? runtime_set_memory
        Int? runtime_set_disk
        Int? runtime_set_max_retries
        Boolean use_ssd = false
    }
    Array[Array[String]]? tsv_to_use = if defined(tsv_of_new_bams) then tsv_of_new_bams else tsv_of_existing_bams
    command <<<
        ## The only purpose of this task is to create the Array[Array[String]] tsv_to_use
    >>>
    output {
        Array[Array[String]] finalized_inputs = tsv_to_use
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
