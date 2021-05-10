version development

task collectScatters {
    input {
        String scatterName
        String scatterIntervals
    }
    command {
    }
    output {
        String allScatterNames = scatterName
        String allScatterIntervals = scatterIntervals
    }
}

task collectAllGroupNames {
    input {
        String groupName
        String sampleName
    }
    command {
    }
    output {
        String allGroupNames = groupName
        String allSampleNames = sampleName
    }
}


task collectExistingBamInfo {
    input {
        String sampleName
        String sampleGroup
        File? sampleBAM
        File? sampleBAI
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
}

task coerceFile {
    input {
        File? input_file
    }
    command <<<
        ## The only purpose of this task is to coerce File? into File
    >>>
    output {
        File coerced = input_file
    }
}
