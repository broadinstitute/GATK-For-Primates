version development

task FinalizeInputs {
    input {
        Array[Array[String]]? tsv_of_new_bams
        Array[Array[String]]? tsv_of_existing_bams
        Int? preemptible_tries
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
        memory: 1 + " GB"
        disks: "local-disk " + 1 + if use_ssd then " SSD" else " HDD"
        preemptible: select_first([preemptible_tries, 5])
    }
}
