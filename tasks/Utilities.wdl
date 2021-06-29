version development

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

##########################################################################
## *** TASK: baseRecalibrator ***
##########################################################################
## Generates before table, applies BQSR, then generates after table.
##########################################################################

import "../structs/structs.wdl"


##########################################################################
## *** TASK: generateBamJson ***
##########################################################################
## Generates a JSON comprising the sample and bam files needed downstream.
##########################################################################

task generateBamJSON {
    input {
        Array[String]? sampleGroups
        Array[String]? sampleNames
        Array[String]? bams
        Array[String]? bam_indexes
        # Runtime
        String container
        Int? runtime_set_preemptible_tries
        Int? runtime_set_cpu
        Int? runtime_set_memory
        Int? runtime_set_disk
        Int? runtime_set_max_retries
        Boolean use_ssd = false
    }
    command <<<
    set -oe pipefail

    python << CODE
    import json

    sampleGroup = ['~{sep="','" sampleGroups}']
    sampleName = ['~{sep="','" sampleNames}']
    bam = ['~{sep="','" bams}']
    bam_index = ['~{sep="','" bam_indexes}']

    data = {}
    data = []

    for i in range(len(sampleGroup)):
        data.append({
            'sampleGroup': sampleGroup[i],
            'sampleName': sampleName[i],
            'bam': bam[i],
            'bam_index': bam_index[i]
        })

    with open('bamJSON.txt', 'w') as outfile:
        json.dump(data, outfile, sort_keys=True, indent=4)

    CODE

    >>>
    output {
        File file = "bamJSON.txt"
    }
    runtime {
        docker: container
        cpu: select_first([runtime_set_cpu, 1])
        gpu: false
        memory: select_first([runtime_set_memory, 1]) + " GB"
        disks: "local-disk " + select_first([runtime_set_disk, 10]) + if use_ssd then " SSD" else " HDD"
        maxRetries: select_first([runtime_set_max_retries, 0])
        preemptible: select_first([runtime_set_preemptible_tries, 5])
        returnCodes: 0
     }
}


##########################################################################
## *** TASK: generateHaplotypeJSON ***
##########################################################################
## Generates a JSON comprising the outputs from HaplotypeCaller, and adds a
## 'Cohort' group containing all samples if the Boolean
## 'multiple_taxonomic_groups' (from ValidateUserInputs.wdl) is true.
##########################################################################

task generateHaplotypeJSON {
    input {
        Array[String] sampleGroups
        Array[String] title_gvcfs
        Array[String] gvcfs
        Array[String] gvcf_indexes
        Boolean multiple_taxonomic_groups
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
    set -oe pipefail

    python << CODE
    import json

    sampleGroup = ['~{sep="','" sampleGroups}']
    title_gvcf = ['~{sep="','" title_gvcfs}']
    gvcf = ['~{sep="','" gvcfs}']
    gvcf_index = ['~{sep="','" gvcf_indexes}']
    make_cohort = '~{multiple_taxonomic_groups}'

    data = {}
    data = []

    if make_cohort == "true":
        for i in range(len(sampleGroup)):
            data.append({
                'groupName': sampleGroup[i],
                'title_gvcf': title_gvcf[i],
                'gvcf': gvcf[i],
                'gvcf_index': gvcf_index[i]
            })
            data.append({
                'groupName': 'Cohort',
                'title_gvcf': title_gvcf[i],
                'gvcf': gvcf[i],
                'gvcf_index': gvcf_index[i]
            })
    else:
        for i in range(len(sampleGroup)):
            data.append({
                'groupName': sampleGroup[i],
                'title_gvcf': title_gvcf[i],
                'gvcf': gvcf[i],
                'gvcf_index': gvcf_index[i]
            })


    with open('haplotypeJSON.txt', 'w') as outfile:
        json.dump(data, outfile, sort_keys=True, indent=4)

    CODE

    >>>
    output {
        File file = "haplotypeJSON.txt"
    }
    runtime {
        docker: container
        cpu: select_first([runtime_set_cpu, 1])
        gpu: false
        memory: select_first([runtime_set_memory, 1]) + " GB"
        disks: "local-disk " + select_first([runtime_set_disk, 10]) + if use_ssd then " SSD" else " HDD"
        maxRetries: select_first([runtime_set_max_retries, 0])
        preemptible: select_first([runtime_set_preemptible_tries, 5])
        returnCodes: 0
     }
}

##########################################################################
## *** TASK: generateGenotypeJSON ***
##########################################################################
## Generates a JSON comprising the genotype data needed downstream.
##########################################################################

task generateGenotypeJSON {
    input {
        Array[String] groupName
        Array[String] loc_unfiltered_vcf
        Array[String] loc_unfiltered_vcf_index
        Array[String] loc_unfiltered_sites_only_vcf
        Array[String] loc_unfiltered_sites_only_vcf_index
        Array[String] loc_hard_filtered_SNPs_vcf
        Array[String] loc_hard_filtered_SNPs_vcf_index
        Array[String] loc_hard_filtered_INDELs_vcf
        Array[String] loc_hard_filtered_INDELs_vcf_index
        Array[String] loc_hard_filtered_SNPs_sites_only_vcf
        Array[String] loc_hard_filtered_SNPs_sites_only_vcf_index
        Array[String] loc_hard_filtered_INDELs_sites_only_vcf
        Array[String] loc_hard_filtered_INDELs_sites_only_vcf_index
        # Runtime
        String container
        Int? runtime_set_preemptible_tries
        Int? runtime_set_cpu
        Int? runtime_set_memory
        Int? runtime_set_disk
        Int? runtime_set_max_retries
        Boolean use_ssd = false
    }
    command <<<
    set -oe pipefail

    python << CODE
    import json

    groupName = ['~{sep="','" groupName}']
    unfiltered_vcf = ['~{sep="','" loc_unfiltered_vcf}']
    unfiltered_vcf_index = ['~{sep="','" loc_unfiltered_vcf_index}']
    unfiltered_sites_only_vcf = ['~{sep="','" loc_unfiltered_sites_only_vcf}']
    unfiltered_sites_only_vcf_index = ['~{sep="','" loc_unfiltered_sites_only_vcf_index}']
    hard_filtered_SNPs_vcf = ['~{sep="','" loc_hard_filtered_SNPs_vcf}']
    hard_filtered_SNPs_vcf_index = ['~{sep="','" loc_hard_filtered_SNPs_vcf_index}']
    hard_filtered_INDELs_vcf = ['~{sep="','" loc_hard_filtered_INDELs_vcf}']
    hard_filtered_INDELs_vcf_index = ['~{sep="','" loc_hard_filtered_INDELs_vcf_index}']
    hard_filtered_SNPs_sites_only_vcf = ['~{sep="','" loc_hard_filtered_SNPs_sites_only_vcf}']
    hard_filtered_SNPs_sites_only_vcf_index = ['~{sep="','" loc_hard_filtered_SNPs_sites_only_vcf_index}']
    hard_filtered_INDELs_sites_only_vcf = ['~{sep="','" loc_hard_filtered_INDELs_sites_only_vcf}']
    hard_filtered_INDELs_sites_only_vcf_index = ['~{sep="','" loc_hard_filtered_INDELs_sites_only_vcf_index}']

    data = {}
    data = []

    for i in range(len(groupName)):
        data.append({
        'groupName': groupName[i],
        'unfiltered_vcf': unfiltered_vcf[i],
        'unfiltered_vcf_index': unfiltered_vcf_index[i],
        'unfiltered_sites_only_vcf': unfiltered_sites_only_vcf[i],
        'unfiltered_sites_only_vcf_index': unfiltered_sites_only_vcf_index[i],
        'hard_filtered_SNPs_vcf': hard_filtered_SNPs_vcf[i],
        'hard_filtered_SNPs_vcf_index': hard_filtered_SNPs_vcf_index[i],
        'hard_filtered_INDELs_vcf': hard_filtered_INDELs_vcf[i],
        'hard_filtered_INDELs_vcf_index': hard_filtered_INDELs_vcf_index[i],
        'hard_filtered_SNPs_sites_only_vcf': hard_filtered_SNPs_sites_only_vcf[i],
        'hard_filtered_SNPs_sites_only_vcf_index': hard_filtered_SNPs_sites_only_vcf_index[i],
        'hard_filtered_INDELs_sites_only_vcf': hard_filtered_INDELs_sites_only_vcf[i],
        'hard_filtered_INDELs_sites_only_vcf_index': hard_filtered_INDELs_sites_only_vcf_index[i]
        })

    with open('genotypeJSON.txt', 'w') as outfile:
        json.dump(data, outfile, sort_keys=True, indent=4)

    CODE

    >>>
    output {
        File file = "genotypeJSON.txt"
    }
    runtime {
        docker: container
        cpu: select_first([runtime_set_cpu, 1])
        gpu: false
        memory: select_first([runtime_set_memory, 1]) + " GB"
        disks: "local-disk " + select_first([runtime_set_disk, 10]) + if use_ssd then " SSD" else " HDD"
        maxRetries: select_first([runtime_set_max_retries, 0])
        preemptible: select_first([runtime_set_preemptible_tries, 5])
        returnCodes: 0
     }
}
