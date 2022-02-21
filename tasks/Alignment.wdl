version 1.0

## Copyright Broad Institute and Graham L Banes, 2021-2022
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
## *** TASK: mapFromPairedFASTQ ***
##########################################################################
## Maps non-interleaved, paired FASTQ inputs to BAM.
##########################################################################

task mapFromPairedFASTQ {
    input {
        File ref
        File ref_dict
        File ref_fai
        Array[File?]+ ref_idxs
        String sampleName
        String? RG_ID
        String? RG_SM
        String? RG_LB
        String? RG_PU
        File? R1
        File? R2
        String execute_aligner
        String? RG_CN
        String? RG_DS
        String? RG_DT
        String? RG_FO
        String? RG_KS
        String? RG_PG
        String? RG_PI
        String? RG_PM
        # Runtime options
        String container
        String path_to_gitc
        Int? runtime_set_preemptible_tries
        Int? runtime_set_cpu
        Int? runtime_set_memory
        Int? runtime_set_disk
        Int? runtime_set_max_retries
        Boolean use_ssd = false
    }
    ## Runtime parameters
    Float size_input_files = size(ref, "GB") + size(ref_dict, "GB") + size(ref_fai, "GB") + size(ref_idxs, "GB") + size(R1, "GB") + size(R2, "GB")
    Int runtime_calculated_disk = ceil(size_input_files * 2.5)
    Int command_mem_gb = select_first([runtime_set_memory, 10]) - 1
    ## Task-specific parameters
    String is_RG_CN = if defined(RG_CN) then "\tCN:~{RG_CN}" else ""
    String is_RG_DS = if defined(RG_DS) then "\tDS:~{RG_DS}" else ""
    String is_RG_DT = if defined(RG_DT) then "\tDT:~{RG_DT}" else ""
    String is_RG_FO = if defined(RG_FO) then "\tFO:~{RG_FO}" else ""
    String is_RG_KS = if defined(RG_KS) then "\tKS:~{RG_KS}" else ""
    String is_RG_PG = if defined(RG_PG) then "\tPG:~{RG_PG}" else ""
    String is_RG_PI = if defined(RG_PI) then "\tPI:~{RG_PI}" else ""
    String is_RG_PM = if defined(RG_PM) then "\tPM:~{RG_PM}" else ""
    String optional_read_groups = "~{is_RG_CN}~{is_RG_DS}~{is_RG_DT}~{is_RG_FO}~{is_RG_KS}~{is_RG_PG}~{is_RG_PI}~{is_RG_PM}"
    command <<<
    set -euo pipefail

        ~{path_to_gitc}~{execute_aligner} -R "@RG\tID:~{RG_ID}\tPL:ILLUMINA\tPU:~{RG_PU}\tLB:~{RG_LB}\tSM:~{RG_SM}~{optional_read_groups}" ~{ref} ~{R1} ~{R2} | samtools view -bt ~{ref_fai} -1 - > ~{sampleName}_mapped_unsorted.bam
        samtools sort -n -o ~{sampleName}_mapped.bam ~{sampleName}_mapped_unsorted.bam

    >>>
    output {
        File output_bam_mapped = "~{sampleName}_mapped.bam"
    }
    runtime {
        docker: container
        cpu: select_first([runtime_set_cpu, 1])
        gpu: false
        memory: select_first([runtime_set_memory, 10]) + " GB"
        disks: "local-disk " + select_first([runtime_set_disk, runtime_calculated_disk]) + if use_ssd then " SSD" else " HDD"
        maxRetries: select_first([runtime_set_max_retries, 0])
        preemptible: select_first([runtime_set_preemptible_tries, 5])
        returnCodes: 0
    }
}

##########################################################################
## *** TASK: mapFromUnmappedBAM ***
##########################################################################
## Streams unmapped BAM to FASTQ, maps to BAM and merges.
##########################################################################

task mapFromUnmappedBAM {
    input {
        File ref
        File ref_dict
        File ref_fai
        Array[File?]+ ref_idxs
        String sampleName
        File? unmapped_bam
        String execute_aligner
        # Runtime options
        String container
        String path_to_gitc
        String path_to_gitc_gatk
        Int? runtime_set_preemptible_tries
        Int? runtime_set_cpu
        Int? runtime_set_memory
        Int? runtime_set_disk
        Int? runtime_set_max_retries
        Boolean use_ssd = false
    }
    ## Runtime parameters
    Float size_input_files = size(ref, "GB") + size(ref_dict, "GB") + size(ref_fai, "GB") + size(ref_idxs, "GB") + size(unmapped_bam, "GB")
    Float disk_multiplier = 2.5
    Int runtime_calculated_disk = ceil(size_input_files + (disk_multiplier * unmapped_bam_size) + 20)
    Int command_mem_gb = select_first([runtime_set_memory, 10]) - 1
    command <<<

        # This is done before "set -o pipefail" because "bwa" will have a rc=1 and we don't want to allow rc=1 to succeed
        # because the sed may also fail with that error and that is something we actually want to fail on.
        BWA_VERSION=$(./bwa 2>&1 | grep -e '^Version' | sed 's/Version: //')

        set -o pipefail
        set -e

        ~{path_to_gitc_gatk}gatk SamToFastq --java-options "-Xmx~{command_mem_gb}G" \
        -I ~{unmapped_bam} \
        -F /dev/stdout \
        --INTERLEAVE true \
        --INCLUDE_NON_PF_READS true \
        | \
        ~{path_to_gitc}~{execute_aligner} ~{ref} /dev/stdin - 2> >(tee ~{sampleName}.bwa.stderr.log >&2) \
        | \
        gatk MergeBamAlignment --java-options "-Xmx~{command_mem_gb}G" \
        --VALIDATION_STRINGENCY SILENT \
        --EXPECTED_ORIENTATIONS FR \
        --ATTRIBUTES_TO_RETAIN X0 \
        --ALIGNED_BAM /dev/stdin \
        --UNMAPPED_BAM ~{unmapped_bam} \
        --OUTPUT ~{sampleName}_mapped_and_merged.bam \
        --REFERENCE_SEQUENCE ~{ref} \
        --PAIRED_RUN true \
        --SORT_ORDER unsorted \
        --IS_BISULFITE_SEQUENCE false \
        --ALIGNED_READS_ONLY false \
        --CLIP_ADAPTERS false \
        --MAX_RECORDS_IN_RAM 2000000 \
        --ADD_MATE_CIGAR true \
        --MAX_INSERTIONS_OR_DELETIONS -1 \
        --PRIMARY_ALIGNMENT_STRATEGY MostDistant \
        --PROGRAM_RECORD_ID "bwamem" \
        --PROGRAM_GROUP_VERSION "$BWA_VERSION" \
        --PROGRAM_GROUP_COMMAND_LINE "~{execute_aligner} ~{ref}" \
        --PROGRAM_GROUP_NAME "bwamem" \
        --UNMAPPED_READ_STRATEGY COPY_TO_TAG \
        --ALIGNER_PROPER_PAIR_FLAGS true \
        --UNMAP_CONTAMINANT_READS true

    >>>
    output {
        File output_mapped_and_merged_bam = "~{sampleName}_mapped_and_merged.bam"
        File output_bam_from_ubam_log = "~{sampleName}.bwa.stderr.log"
    }
    runtime {
        docker: container
        cpu: select_first([runtime_set_cpu, 1])
        gpu: false
        memory: select_first([runtime_set_memory, 10]) + " GB"
        disks: "local-disk " + select_first([runtime_set_disk, runtime_calculated_disk]) + if use_ssd then " SSD" else " HDD"
        maxRetries: select_first([runtime_set_max_retries, 0])
        preemptible: select_first([runtime_set_preemptible_tries, 5])
        returnCodes: 0
    }
}
