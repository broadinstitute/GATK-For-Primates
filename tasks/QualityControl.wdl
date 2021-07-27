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

import "../structs/structs.wdl"

##########################################################################
## *** TASK: failWithError ***
##########################################################################
## Fails the workflow with the given error message written to stderr.
##########################################################################

task failWithError {
    input {
        String message = message
    }
    command <<<
    set -euo pipefail

    python <<CODE

    import sys

    message = "~{message}"

    sys.stderr.write(message)
    sys.exit(1)

    CODE
    >>>
}

##########################################################################
## *** TASK: validateVCF ***
##########################################################################
## Validates input VCF file for conformance to VCF standard.
##########################################################################

task validateVCF {
    input {
        File ref
        File ref_dict
        File ref_fai
        File? variant_file
        File? variant_file_index
        String? optional_arguments
        # Runtime options
        String container
        Int? runtime_set_preemptible_tries
        Int? runtime_set_cpu
        Int? runtime_set_memory
        Int? runtime_set_disk
        Int? runtime_set_max_retries
        Boolean use_ssd = false
    }
    ## Runtime parameters
    Float size_input_files = size(ref, "GB") + size(ref_dict, "GB") + size(ref_fai, "GB") + size(variant_file, "GB") + size(variant_file_index, "GB")
    Int runtime_calculated_disk = ceil(size_input_files * 1) + 25
    Int command_mem_gb = select_first([runtime_set_memory, 4]) - 1
    ## Task parameters
    String reference_flag = if defined(ref) then "-R ~{ref} " else ""
    String variant_flag = if defined(variant_file) then "-V ~{variant_file} " else ""
    command <<<
    set -oe pipefail

        gatk ValidateVariants --java-options -Xmx~{command_mem_gb}G \
        -R ~{ref} \
        -V ~{variant_flag} \
        ~{optional_arguments}

    >>>
    runtime {
        docker: container
        cpu: select_first([runtime_set_cpu, 1])
        gpu: false
        memory: select_first([runtime_set_memory, 4]) + " GB"
        disks: "local-disk " + select_first([runtime_set_disk, 3]) + if use_ssd then " SSD" else " HDD"
        maxRetries: select_first([runtime_set_max_retries, 0])
        preemptible: select_first([runtime_set_preemptible_tries, 5])
        returnCodes: 0
     }
}

##########################################################################
## *** TASK: validateRecords ***
##########################################################################
## Checks for errors in input variables provided for each sample.
##########################################################################

task validateRecords {
    input {
        String groupName
        String sampleName
        String? R1
        String? R2
        String? RG_ID
        String? RG_SM
        String? RG_LB
        String? RG_PU
        String? bam
        String? bam_index
        String? unmapped_bam
    }
    command <<<
    set -euo pipefail

    python <<CODE

    import sys
    
    groupName = "~{groupName}"
    sampleName = "~{sampleName}"
    R1 = "~{R1}"
    R2 = "~{R2}"
    RG_ID = "~{RG_ID}"
    RG_SM = "~{RG_SM}"
    RG_LB = "~{RG_LB}"
    RG_PU = "~{RG_PU}"
    bam = "~{bam}"
    bam_index = "~{bam_index}"
    unmapped_bam = "~{unmapped_bam}"

    # If no groupName is provided
    if groupName == "NULL":
        sys.stderr.write("At least one sample is missing a 'taxon_group' input variable.")
        sys.exit(1)
    else:
        # If groupName matches 'cohort'
        if groupName.casefold() == "cohort":
            sys.stderr.write("The taxon_group 'Cohort' is reserved for use by the workflow.")
            sys.exit(1)

    # If no sampleName is provided
    if sampleName == "NULL":
        sys.stderr.write("At least one sample is missing a 'name' input variable.")
        sys.exit(1)

    # If unmapped_bam is provided
    if unmapped_bam != "NULL":
        if R1 != "NULL":
            sys.stderr.write("You must provide either FASTQ reads (i.e. R1 and R2) or an 'unmapped_bam', but not both.")
            sys.exit(1)
        if R2 != "NULL":
            sys.stderr.write("You must provide either FASTQ reads (i.e. R1 and R2) or an 'unmapped_bam', but not both.")
            sys.exit(1)
        if RG_ID != "NULL":
            sys.stderr.write("Read group information (i.e. RG_ID) cannot be provided alongside an 'unmapped_bam'. Ensure all read group information is contained within the bam file.")
            sys.exit(1)
        if RG_SM != "NULL":
            sys.stderr.write("Read group information (i.e. RG_SM) cannot be provided alongside an 'unmapped_bam'. Ensure all read group information is contained within the bam file.")
            sys.exit(1)
        if RG_LB != "NULL":
            sys.stderr.write("Read group information (i.e. RG_LB) cannot be provided alongside an 'unmapped_bam'. Ensure all read group information is contained within the bam file.")
            sys.exit(1)
        if RG_PU != "NULL":
            sys.stderr.write("Read group information (i.e. RG_PU) cannot be provided alongside an 'unmapped_bam'. Ensure all read group information is contained within the bam file.")
            sys.exit(1)
        if bam != "NULL":
            sys.stderr.write("You must provide either a mapped bam (i.e. 'bam') or an 'unmapped_bam', but not both.")
            sys.exit(1)
        if bam_index != "NULL":
            sys.stderr.write("You must provide either a mapped bam index (i.e. 'bam_index') or an 'unmapped_bam', but not both.")
            sys.exit(1)

    # If bam is provided
    if bam != "NULL":
        if bam_index == "NULL":
            sys.stderr.write("You must provide a bam_index for each bam file.")
            sys.exit(1)
        if R1 != "NULL":
            sys.stderr.write("You must provide either FASTQ reads (i.e. R1 and R2) or a mapped 'bam', but not both.")
            sys.exit(1)
        if R2 != "NULL":
            sys.stderr.write("You must provide either FASTQ reads (i.e. R1 and R2) or a mapped 'bam', but not both.")
            sys.exit(1)
        if RG_ID != "NULL":
            sys.stderr.write("Read group information (i.e. RG_ID) cannot be provided alongside a 'bam'. Ensure all read group information is contained within the bam file.")
            sys.exit(1)
        if RG_SM != "NULL":
            sys.stderr.write("Read group information (i.e. RG_SM) cannot be provided alongside a 'bam'. Ensure all read group information is contained within the bam file.")
            sys.exit(1)
        if RG_LB != "NULL":
            sys.stderr.write("Read group information (i.e. RG_LB) cannot be provided alongside a 'bam'. Ensure all read group information is contained within the bam file.")
            sys.exit(1)
        if RG_PU != "NULL":
            sys.stderr.write("Read group information (i.e. RG_PU) cannot be provided alongside a 'bam'. Ensure all read group information is contained within the bam file.")
            sys.exit(1)

    # If only one FASTQ is provided
    if R1 != "NULL":
        if R2 == "NULL":
            sys.stderr.write("You must provide both an R1 and an R2 as paired-end non-interleaved inputs.")
            sys.exit(1)
    if R2 != "NULL":
        if R1 == "NULL":
            sys.stderr.write("You must provide both an R1 and an R2 as paired-end non-interleaved inputs.")
            sys.exit(1)

    # If FASTQ files are provided
    if R1 != "NULL":
        if RG_ID == "NULL":
            sys.stderr.write("Read group information (i.e. RG_ID) must be provided alongside FASTQ reads.")
            sys.exit(1)

        if RG_SM == "NULL":
            sys.stderr.write("Read group information (i.e. RG_ID) must be provided alongside FASTQ reads.")
            sys.exit(1)

        if RG_LB == "NULL":
            sys.stderr.write("Read group information (i.e. RG_ID) must be provided alongside FASTQ reads.")
            sys.exit(1)

        if RG_PU == "NULL":
            sys.stderr.write("Read group information (i.e. RG_ID) must be provided alongside FASTQ reads.")
            sys.exit(1)

    CODE

    >>>
    output {
        String groupNames = "~{groupName}"
        String sampleNames = "~{sampleName}"
        String R1s = "~{R1}"
        String R2s = "~{R2}"
        String RG_IDs = "~{RG_ID}"
        String RG_SMs = "~{RG_SM}"
        String RG_LBs = "~{RG_LB}"
        String RG_PUs = "~{RG_PU}"
        String bams = "~{bam}"
        String bam_indexes = "~{bam_index}"
        String unmapped_bams = "~{unmapped_bam}"
    }
}

##########################################################################
## *** TASK: validateCohort ***
##########################################################################
## Checks for duplicate input variables and multiple taxa.
##########################################################################

task validateCohort {
    input {
        Array[String] scatterNames
        Array[String] scatterIntervals
        Array[String] groupNames
        Array[String] sampleNames
        Array[String] RG_IDs
        Array[String] RG_SMs
        Array[String] RG_LBs
        # Runtime options
        String container
        Int? runtime_set_preemptible_tries
        Int? runtime_set_cpu
        Int? runtime_set_memory
        Int? runtime_set_disk
        Int? runtime_set_max_retries
        Boolean use_ssd = false
    }
    command <<<
    set -euo pipefail

    python <<CODE

    import sys

    # Import WDL arrays to Python
    scatterName = ['~{sep="','" scatterNames}']
    scatterInterval = ['~{sep="','" scatterIntervals}']
    groupName = ['~{sep="','" groupNames}']
    sampleName = ['~{sep="','" sampleNames}']
    RG_ID = ['~{sep="','" RG_IDs}']
    RG_SM = ['~{sep="','" RG_SMs}']
    RG_LB = ['~{sep="','" RG_LBs}']

    # Remove NULL values from read groups
    RG_ID[:] = (value for value in RG_ID if value != "NULL")
    RG_SM[:] = (value for value in RG_SM if value != "NULL")
    RG_LB[:] = (value for value in RG_LB if value != "NULL")

    # Get unique values only
    unique_scatterNames = list(dict.fromkeys(scatterName))
    unique_scatterIntervals = list(dict.fromkeys(scatterInterval))
    unique_groupNames = list(dict.fromkeys(groupName))
    unique_sampleNames = list(dict.fromkeys(sampleName))
    unique_RG_IDs = list(dict.fromkeys(RG_ID))
    unique_RG_SMs = list(dict.fromkeys(RG_SM))
    unique_RG_LBs = list(dict.fromkeys(RG_LB))

    # Determine if multiple taxa are being used
    if len(groupName) is not len(unique_groupNames):
        sys.stdout.write("true")
    else:
        sys.stdout.write("false")

    # Check for duplicate scatter names
    if len(scatterName) is not len(unique_scatterNames):
        sys.stderr.write("You appear to have duplicate scatter names. Every scatter must have a unique name.")
        sys.exit(1)

    # Check that every scatter has an interval
    if len(scatterInterval) is not len(unique_scatterIntervals):
        sys.stderr.write("You appear to have scatters that are missing either a name or a list of intervals. Each defined scatter must have both.")
        sys.exit(1)

    # Check for duplicate sample names
    if len(sampleName) is not len(unique_sampleNames):
        sys.stderr.write("You appear to have duplicate sample names. If you have multiple FASTQ files from the same sample, you must submit these data in unmapped BAM format.")
        sys.exit(1)

    # Check for duplicate RG_IDs
    if len(RG_ID) is not len(unique_RG_IDs):
        sys.stderr.write("You appear to have duplicate RG_IDs.")
        sys.exit(1)

    # Check for duplicate RG_SMs
    if len(RG_SM) is not len(unique_RG_SMs):
        sys.stderr.write("You appear to have duplicate RG_SMs. If you have multiple FASTQ files from the same sample, you must submit these data in unmapped BAM format.")
        sys.exit(1)

    # Check for duplicate RG_LBs
    if len(RG_LB) is not len(unique_RG_LBs):
        sys.stderr.write("You appear to have duplicate RG_LBs. If you have multiple FASTQ files from the same library, you must submit these data in unmapped BAM format.")
        sys.exit(1)

    CODE

    >>>
    output {
        Array[String] output_scatterNames = scatterNames
        Boolean multiple_taxonomic_groups = read_boolean(stdout())
    }
}
