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
## *** TASK: unpackagePolymorphicRegions ***
##########################################################################
## Unpackages `packaged_polymorphic_regions` in 'repeat' and 'final'
## modes, generates blank JSON in 'initial' mode to avoid Cromwell failure
##########################################################################

task unpackagePolymorphicRegions {
    input {
        String mode
        File? packaged_polymorphic_regions
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
    set -euo pipefail

    python <<CODE

    import sys
    import json
    import tarfile

    mode = "~{mode}"

    if mode == "initial":
        data = {}
        data = []
        data.append({
            'scatterName': 'NULL',
            'intervalList': 'NULL'
        })
        with open('polymorphicRegions.json', 'w') as outfile:
            json.dump(data, outfile, sort_keys=True, indent=4)

    if mode != "initial":
        try:
            packaged_tar = tarfile.open('~{packaged_polymorphic_regions}')
            packaged_tar.extractall('.')
            packaged_tar.close()
        except FileNotFoundError:
            sys.stderr.write("Cannot locate or read packaged polymorphic regions.")
            sys.exit(1)

    CODE
    >>>
    output {
        File polymorphicRegionsJSON = "polymorphicRegions.json"
        Array[File] polymorphicRegionsIntervalLists = glob("*.interval_list")
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
## *** TASK: packagePolymorphicRegions ***
##########################################################################
## Generates a JSON comprising the filenames of the interval lists per
## scatter, and packages this into a .tar.gz with the interval lists.
##########################################################################

task packagePolymorphicRegions {
    input {
        Array[String] scatterNames
        Array[File] intervalLists
        Array[String] intervalListFilenames
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
    import tarfile
    import os.path

    scatterName = ['~{sep="','" scatterNames}']
    intervalList = ['~{sep="','" intervalLists}']
    intervalListFilenames = ['~{sep="','" intervalListFilenames}']

    data = {}
    data = []

    for i in range(len(scatterName)):
        data.append({
            'scatterName': scatterName[i],
            'intervalList': intervalListFilenames[i],
        })

    with open('polymorphicRegions.json', 'w') as outfile:
        json.dump(data, outfile, sort_keys=True, indent=4)

    tar = tarfile.open("polymorphicRegions.tar.gz", "w:gz")
    for name in ['~{sep="','" intervalLists}']:
        tar.add(name, arcname=os.path.basename(name))
    tar.add('polymorphicRegions.json')
    tar.close()

    CODE

    >>>
    output {
        File file = "polymorphicRegions.json"
        File file_tar = "polymorphicRegions.tar.gz"
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
        Array[String] gvcfs
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
    gvcf = ['~{sep="','" gvcfs}']
    make_cohort = '~{multiple_taxonomic_groups}'

    data = {}
    data = []

    if make_cohort == "true":
        for i in range(len(sampleGroup)):
            data.append({
                'groupName': sampleGroup[i],
                'gvcf': gvcf[i]
            })
            data.append({
                'groupName': 'Cohort',
                'gvcf': gvcf[i]
            })
    else:
        for i in range(len(sampleGroup)):
            data.append({
                'groupName': sampleGroup[i],
                'gvcf': gvcf[i]
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

##############################################################################################################################################################################################################################

##########################################################################
## *** TASK: pipeHaplotypes ***
##########################################################################
## Searches a JSON struct for the query.
##########################################################################

task pipeHaplotypes {
    input {
        String groupName
        String scatterName
        File haplotypeJSON
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

    with open('~{haplotypeJSON}') as f:
        haplotypeJSON = json.load(f)

    with open('gvcfs.txt', 'w') as f_gvcf, open('gvcf_indexes.txt', 'w') as f_gvcf_indexes:
        for keyval in haplotypeJSON:
            if keyval['groupName'] == '~{groupName}' and keyval['scatterName'] == '~{scatterName}':
                f_gvcf.write('keyval['gvcf']' + '\n')
                f_gvcf_indexes.write('keyval['gvcf_index']' + '.tbi\n')

    CODE

    >>>
    output {
        Array[String] gvcfs = read_lines("gvcfs.txt")
        Array[String] gvcf_indexes = read_lines("gvcf_indexes.txt")
    }
    runtime {
        docker: container
        cpu: select_first([runtime_set_cpu, 1])
        gpu: false
        memory: select_first([runtime_set_memory, 2]) + " GB"
        disks: "local-disk " + select_first([runtime_set_disk, 5]) + if use_ssd then " SSD" else " HDD"
        maxRetries: select_first([runtime_set_max_retries, 0])
        preemptible: select_first([runtime_set_preemptible_tries, 5])
        returnCodes: 0
     }

}




##########################################################################
## *** TASK: makeBlankJSON ***
##########################################################################
## Makes a blank JSON file.
##########################################################################

task makeBlankJSON {
    input {
        # Runtime
        String container
    }
    command <<<
    set -oe pipefail

    python << CODE
    import json

    data = {}
    data = []

    data.append({
        'scatterName': 'NULL',
        'intervalList': 'NULL'
    })

    with open('blankJSON.txt', 'w') as outfile:
        json.dump(data, outfile, sort_keys=True, indent=4)

    CODE

    >>>
    output {
        File file = "blankJSON.txt"
    }
    runtime {
        docker: container
        cpu: "1"
        gpu: false
        memory: "1 GB"
        disks: "local-disk " + "5 HDD"
        maxRetries: "0"
        preemptible: "5"
        returnCodes: 0
     }
}






##########################################################################
## *** TASK: fetchIntervalList ***
##########################################################################
## Makes a blank JSON file.
##########################################################################

task fetchIntervalListByScatter {
    input {
        # Runtime
        String container
    }
    command <<<
    set -oe pipefail

    python << CODE
    import json

    data = {}
    data = []

    data.append({
        'scatterName': 'NULL',
        'intervalList': 'NULL'
    })

    with open('blankJSON.txt', 'w') as outfile:
        json.dump(data, outfile, sort_keys=True, indent=4)

    CODE

    >>>
    output {
        File file = "blankJSON.txt"
    }
    runtime {
        docker: container
        cpu: "1"
        gpu: false
        memory: "1 GB"
        disks: "local-disk " + "5 HDD"
        maxRetries: "0"
        preemptible: "5"
        returnCodes: 0
     }
}
