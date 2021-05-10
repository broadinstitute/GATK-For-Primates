version development

task FinalizeInputs {
    input {
        Array[Array[String]]? tsv_of_new_bams
        Array[Array[String]]? tsv_of_existing_bams
    }
    Array[Array[String]]? map_to_use = if defined(tsv_of_new_bams) then tsv_of_new_bams else tsv_of_existing_bams
    command <<<
        ## The only purpose of this task is to create the Array[Array[String]] map_to_use
    >>>
    output {
        Array[Array[String]] finalized_inputs = map_to_use
    }
}