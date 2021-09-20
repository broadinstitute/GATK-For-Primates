version 1.0

## Copyright Broad Institute and Wisconsin National Primate Research Center,
## University of Wisconsin-Madison, 2021
## 
##########################################################################
## TAKE NOTE: This WDL acts as a wrapper for the 'GATK for Primates'
## pipeline to run on Terra (https://terra.bio). You should not use this
## WDL to run the pipeline anywhere else: it will fail! Instead, use:
## gatk-for-primates-germline-snps-indels.wdl
##########################################################################
##
## Complete germline short variant discovery pipeline optimized for
## non-human primates, following proposed GATK Best Practices for
## Non-Human Animal Genomes.
##
## For requirements, expectations and outputs, please review the complete
## documentation at: https://github.com/broadinstitute/GATK-For-Primates
##
## Software version requirements :
## - Cromwell 67
## - bwa 0.7.15 (note: from GITC)
## - Samtools 1.11 (note: from GITC)
## - GATK 4.2.1.0 (note: GATK 4.1.8.0 is used in GITC)
## - Python 3.9.5
##
## Program versions can be changed by defining alternative containers.
## Runtime parameters are optimized for Terra (https://terra.bio/)
##
## LICENSING : 
## This script is released under the WDL source code license (BSD-3) (see LICENSE in 
## https://github.com/broadinstitute/wdl). However, the programs it calls may be
## subject to different licenses. Users are responsible for checking that they are
## authorized to run all programs before running this script. Please see the docker
## containers for detailed licensing information pertaining to the included programs.

import "gatk-for-primates-germline-snps-indels.wdl" as gatkForPrimates
import "structs/structs.wdl"

##########################################################################
## WORKFLOW DEFINITION
##########################################################################

workflow GATKForPrimatesOnTerra {

    ##########################################################################
    ## DEFINE WORFKLOW INPUTS
    ##########################################################################

    input {

        ## Collect workflow mode from Terra
        String mode # options: initial / repeat / final
 
        ## Collect optional variables from Terra
        Boolean validate_truth_sets = true # options: true / false; if false this will disable running ValidateVariants on truth sets in 'Final' mode
        Boolean flowcell_patterned = true # options: true / false; this influences pixel distance when marking duplicates
        Int? merge_contigs_into_num_partitions # options: optional parameter for GenomicsDBImport
        Boolean bwamem2 = false # options: true / false; indicating bwa (as bwa mem) or bwamem2 (as bwamem2 mem) ***-Coming-Soon-***
        Boolean cram_not_bam = true # options: true / false; if false this will disable use of CRAM instead of BAM format

        ## Collect Terra workspace and project information
        String workspace_name
        String terra_project

        ## Collect reference files from Terra
        File ref
        File ref_dict
        File ref_fai
        File? ref_amb
        File? ref_ann
        File? ref_pac
        File? ref_bwt
        File? ref_sa
        File? ref_0123
        File? ref_bwt_2bit_64

        ## Collect scatterList JSON file from Terra
        File scatterList_json

        ## Collect .tar.gz of packaged interval lists from 'initial' mode
        File? packaged_polymorphic_regions

        ## Collect truth sets from Terra
        File? truth_set_SNPs # options: optional SNP truth set (sites-only VCF file) for VQSR (training set is produced via hard filtering)
        File? truth_set_SNPs_index # index for the above
        File? truth_set_INDELs # options: optional INDEL truth set (sites-only VCF file) for VQSR (training set is produced via hard filtering)
        File? truth_set_INDELs_index # index for the above
        
        ## Define arrays from input JSONs; definitions are in the structs/structs.wdl file
        Array[scatterInfo]+ scatterList = read_json(scatterList_json)
        
        ## Define containers
        String container_gatk = "broadinstitute/gatk:4.2.1.0"
        String container_gitc = "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.5.7-2021-06-09_16-47-48Z"
        String container_python = "python:3.9.5"

        ## Define paths
        String path_to_gitc = "/usr/gitc/"
        String path_to_gitc_gatk = "/usr/gitc/gatk4/"

        ## TERRA ONLY
        Array[String] name
        Array[String]? RG_ID
        Array[String]? RG_LB
        Array[String]? RG_SM
        Array[String]? RG_PU
        Array[String] taxon_group
        Array[String]? R1
        Array[String]? R2
        Array[String]? unmapped_bam
        Array[String]? bam
        Array[String]? bam_index
        Array[String]? RG_CN
        Array[String]? RG_DS
        Array[String]? RG_DT
        Array[String]? RG_FO
        Array[String]? RG_KS
        Array[String]? RG_PG
        Array[String]? RG_PI
        Array[String]? RG_PM        

    }

    call generateSampleJSONforTerra as generateSampleJSONforTerra {
        input:
            name = name,
            RG_ID = RG_ID,
            RG_LB = RG_LB,
            RG_SM = RG_SM,
            RG_PU = RG_PU,
            taxon_group = taxon_group,
            R1 = R1,
            R2 = R2,
            unmapped_bam = unmapped_bam,
            bam = bam,
            bam_index = bam_index,
            RG_CN = RG_CN,
            RG_DS = RG_DS,
            RG_DT = RG_DT,
            RG_FO = RG_FO,
            RG_KS = RG_KS,
            RG_PG = RG_PG,
            RG_PI = RG_PI,
            RG_PM = RG_PM, 
            container = container_python,
    }

    Array[sampleInfo]+ sampleList = read_json(generateSampleJSONforTerra.file)

    call gatkForPrimates.GATKForPrimatesGermlineSNPsIndels_GATK4 as gatkForPrimates {
        input:
            mode = mode,
            validate_truth_sets = validate_truth_sets,
            flowcell_patterned = flowcell_patterned,
            merge_contigs_into_num_partitions = merge_contigs_into_num_partitions,
            bwamem2 = bwamem2,
            cram_not_bam = cram_not_bam,
            ref = ref,
            ref_dict = ref_dict,
            ref_fai = ref_fai,
            ref_amb = ref_amb,
            ref_ann = ref_ann,
            ref_pac = ref_pac,
            ref_bwt = ref_bwt,
            ref_sa = ref_sa,
            ref_0123 = ref_0123,
            ref_bwt_2bit_64 = ref_bwt_2bit_64,
            packaged_polymorphic_regions = packaged_polymorphic_regions,
            truth_set_SNPs = truth_set_SNPs,
            truth_set_SNPs_index = truth_set_SNPs_index,
            truth_set_INDELs = truth_set_INDELs,
            truth_set_INDELs_index = truth_set_INDELs_index,
            scatterList = scatterList,
            sampleList = sampleList,
            container_gatk = container_gatk,
            container_gitc = container_gitc,
            container_python = container_python,
            path_to_gitc = path_to_gitc,
            path_to_gitc_gatk = path_to_gitc_gatk,
    }
    
    if (mode != "final") {
        call collectTerraOutputs {
            input:
                recalibrated_bam = gatkForPrimates.recalibrated_bam,
                recalibrated_bam_index = gatkForPrimates.recalibrated_bam_index,
                table_before = gatkForPrimates.table_before,
                table_after = gatkForPrimates.table_after,
                plots = gatkForPrimates.plots,
                recalibrated_sampleName = gatkForPrimates.recalibrated_sampleName,
                container = container_python,
        }
        call upsertToTerra {
            input:
                tsv_file = collectTerraOutputs.tsv_to_upsert,
                workspace_name = workspace_name,
                terra_project = terra_project,
                container = "schaluvadi/pathogen-genomic-surveillance:batch_upsert",
        }
    }
    
    output {
        ## Outputs from Terra wrapper:
        File? terraJSON = generateSampleJSONforTerra.file
        File? tsv_to_upsert = collectTerraOutputs.tsv_to_upsert
        File? upsert_json = upsertToTerra.upsert_json

        ## Outputs from initial mode:
        File? polymorphic_sites_JSON = gatkForPrimates.polymorphic_sites_JSON
        File? polymorphic_sites_tar = gatkForPrimates.polymorphic_sites_tar
        ## Outputs from initial/repeat modes:
        File? high_confidence_sites_BQSR_vcf = gatkForPrimates.high_confidence_sites_BQSR_vcf
        File? high_confidence_sites_BQSR_vcf_index = gatkForPrimates.high_confidence_sites_BQSR_vcf_index
        Array[File]? recalibrated_bam = gatkForPrimates.recalibrated_bam
        Array[File]? recalibrated_bam_index = gatkForPrimates.recalibrated_bam_index
        Array[File]? table_before = gatkForPrimates.table_before
        Array[File]? table_after = gatkForPrimates.table_after
        Array[File]? plots = gatkForPrimates.plots
        ## Outputs from final mode
        ##File? finalGenotypes = mergeFinalGenotypes.output_merged_vcf
        ##File? finalGenotypesIndex = mergeFinalGenotypes.output_merged_vcf_index
        ##File? finalCallset = gatherFinalCallset.output_vcf
        ##File? finalCallsetIndex = gatherFinalCallset.output_vcf_index

    }
    
}

##########################################################################
## *** TASK: generateSampleJSONforTerra ***
##########################################################################
## Generates a workflow-required sample JSON from the Terra inputs.
##########################################################################

task generateSampleJSONforTerra {
    input {
        Array[String] name
        Array[String]? RG_ID
        Array[String]? RG_LB
        Array[String]? RG_SM
        Array[String]? RG_PU
        Array[String] taxon_group
        Array[String]? R1
        Array[String]? R2
        Array[String]? unmapped_bam
        Array[String]? bam
        Array[String]? bam_index
        Array[String]? RG_CN
        Array[String]? RG_DS
        Array[String]? RG_DT
        Array[String]? RG_FO
        Array[String]? RG_KS
        Array[String]? RG_PG
        Array[String]? RG_PI
        Array[String]? RG_PM
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

    name = ['~{sep="','" name}']
    RG_ID = ['~{sep="','" RG_ID}']
    RG_LB = ['~{sep="','" RG_LB}']
    RG_SM = ['~{sep="','" RG_SM}']
    RG_PU = ['~{sep="','" RG_PU}']
    taxon_group = ['~{sep="','" taxon_group}']
    R1 = ['~{sep="','" R1}']
    R2 = ['~{sep="','" R2}']
    unmapped_bam = ['~{sep="','" unmapped_bam}']
    bam = ['~{sep="','" bam}']
    bam_index = ['~{sep="','" bam_index}']
    RG_CN = ['~{sep="','" RG_CN}']
    RG_DS = ['~{sep="','" RG_DS}']
    RG_DT = ['~{sep="','" RG_DT}']
    RG_FO = ['~{sep="','" RG_FO}']
    RG_KS = ['~{sep="','" RG_KS}']
    RG_PG = ['~{sep="','" RG_PG}']
    RG_PI = ['~{sep="','" RG_PI}']
    RG_PM = ['~{sep="','" RG_PM}']

    data = {}
    data = []

    for i in range(len(name)):

        try:
            out_R1 = R1[i]
        except IndexError:
            out_R1 = "NULL"

        try:
            out_R1 = R1[i]
        except IndexError:
            out_R1 = "NULL"

        try:
            out_R2 = R2[i]
        except IndexError:
            out_R2 = "NULL"

        try:
            out_unmapped_bam = unmapped_bam[i]
        except IndexError:
            out_unmapped_bam = "NULL"

        try:
            out_bam = bam[i]
        except IndexError:
            out_bam = "NULL"

        try:
            out_bam_index = bam_index[i]
        except IndexError:
            out_bam_index = "NULL"

        try:
            out_RG_CN = RG_CN[i]
        except IndexError:
            out_RG_CN = "NULL"

        try:
            out_RG_DS = RG_DS[i]
        except IndexError:
            out_RG_DS = "NULL"

        try:
            out_RG_DT = RG_DT[i]
        except IndexError:
            out_RG_DT = "NULL"

        try:
            out_RG_FO = RG_FO[i]
        except IndexError:
            out_RG_FO = "NULL"

        try:
            out_RG_KS = RG_KS[i]
        except IndexError:
            out_RG_KS = "NULL"

        try:
            out_RG_PG = RG_PG[i]
        except IndexError:
            out_RG_PG = "NULL"

        try:
            out_RG_PI = RG_PI[i]
        except IndexError:
            out_RG_PI = "NULL"

        try:
            out_RG_PM = RG_PM[i]
        except IndexError:
            out_RG_PM = "NULL"

        data.append({
            'name': name[i],
            'RG_ID': RG_ID[i],
            'RG_LB': RG_LB[i],
            'RG_SM': RG_SM[i],
            'RG_PU': RG_PU[i],
            'taxon_group': taxon_group[i],
            'R1': out_R1,
            'R2': out_R2,
            'unmapped_bam': out_unmapped_bam,
            'bam': out_bam,
            'bam_index': out_bam_index,
            'RG_CN': out_RG_CN,
            'RG_DS': out_RG_DS,
            'RG_DT': out_RG_DT,
            'RG_FO': out_RG_FO,
            'RG_KS': out_RG_KS,
            'RG_PG': out_RG_PG,
            'RG_PI': out_RG_PI,
            'RG_PM': out_RG_PM,
            })

    with open('terraInputs_prep.txt', 'w') as tempfile:
        json.dump(data, tempfile, sort_keys=True, indent=4)

    obj = json.load(open("terraInputs_prep.txt"))

    for n in range(len(obj)):
        print(n)

    for i in range(len(obj)):
        if obj[i]["R1"] == "NULL" or obj[i]["R1"] == "":
            obj[i].pop("R1")
        if obj[i]["R2"] == "NULL" or obj[i]["R2"] == "":
            obj[i].pop("R2")
        if obj[i]["unmapped_bam"] == "NULL" or obj[i]["unmapped_bam"] == "":
            obj[i].pop("unmapped_bam")
        if obj[i]["bam"] == "NULL" or obj[i]["bam"] == "":
            obj[i].pop("bam")
        if obj[i]["bam_index"] == "NULL" or obj[i]["bam_index"] == "":
            obj[i].pop("bam_index")
        if obj[i]["RG_CN"] == "NULL" or obj[i]["RG_CN"] == "":
            obj[i].pop("RG_CN")
        if obj[i]["RG_DS"] == "NULL" or obj[i]["RG_DS"] == "":
            obj[i].pop("RG_DS")
        if obj[i]["RG_DT"] == "NULL" or obj[i]["RG_DT"] == "":
            obj[i].pop("RG_DT")
        if obj[i]["RG_FO"] == "NULL" or obj[i]["RG_FO"] == "":
            obj[i].pop("RG_FO")
        if obj[i]["RG_KS"] == "NULL" or obj[i]["RG_KS"] == "":
            obj[i].pop("RG_KS")
        if obj[i]["RG_PG"] == "NULL" or obj[i]["RG_PG"] == "":
            obj[i].pop("RG_PG")
        if obj[i]["RG_PI"] == "NULL" or obj[i]["RG_PI"] == "":
            obj[i].pop("RG_PI")
        if obj[i]["RG_PM"] == "NULL" or obj[i]["RG_PM"] == "":
            obj[i].pop("RG_PM")    

    with open('terraInputs.txt', 'w') as outfile:
        json.dump(obj, outfile, sort_keys=True, indent=4)

    CODE

    >>>
    output {
        File file = "terraInputs.txt"
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
## *** TASK: collectTerraOutputs ***
##########################################################################
## Creates a TSV of sample outputs.
##########################################################################

task collectTerraOutputs {
    input {
        Array[String]? recalibrated_bam
        Array[String]? recalibrated_bam_index
        Array[String]? table_before
        Array[String]? table_after
        Array[String]? plots
        Array[String]? recalibrated_sampleName
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
    recalibrated_bam = ['~{sep="','" recalibrated_bam}']
    recalibrated_bam_index = ['~{sep="','" recalibrated_bam_index}']
    table_before = ['~{sep="','" table_before}']
    table_after = ['~{sep="','" table_after}']
    plots = ['~{sep="','" plots}']
    recalibrated_sampleName = ['~{sep="','" recalibrated_sampleName}']

    if len(recalibrated_bam)!= len(recalibrated_bam_index) != len(table_before) != len(table_after) != len(plots) != len(recalibrated_sampleName):
        print("Numbers of input variables are not equal.")
        exit(1)

    with open("tsv_to_upsert.tsv", "w") as fi:
        fi.write("entity:sample_id\trecalibrated_bam\trecalibrated_bam_index\ttable_before\ttable_after\tplots\n")
        for i in range(len(recalibrated_bam)):
            fi.write(recalibrated_sampleName[i] + "\t" + recalibrated_bam[i] + "\t" + recalibrated_bam_index[i] + "\t" + table_before[i] + "\t" + table_after[i] + "\t" + plots[i] + "\n")
    CODE
    >>>
    runtime {
        docker: container
        cpu: select_first([runtime_set_cpu, 1])
        gpu: false
        memory: select_first([runtime_set_memory, 2]) + " GB"
        disks: "local-disk " + select_first([runtime_set_disk, 10]) + if use_ssd then " SSD" else " HDD"
        maxRetries: select_first([runtime_set_max_retries, 0])
        preemptible: select_first([runtime_set_preemptible_tries, 5])
        returnCodes: 0
     }
    output {
        File tsv_to_upsert = "tsv_to_upsert.tsv"
    }
}

##########################################################################
## *** TASK: upsertToTerra ***
##########################################################################
## 'Upserts' entities to Terra.
##########################################################################

task upsertToTerra {
    input {
        File tsv_file
        String workspace_name
        String terra_project
        # Runtime
        String container
        Int? runtime_set_preemptible_tries
        Int? runtime_set_cpu
        Int? runtime_set_memory
        Int? runtime_set_disk
        Int? runtime_set_max_retries
        Boolean use_ssd = false
    }
    command {
    set -e
    wget https://raw.githubusercontent.com/broadinstitute/GATK-For-Primates/main/scripts/upsert_entities.py
    python3 upsert_entities.py \
        -t "~{tsv_file}" \
        -p "~{terra_project}" \
        -w "~{workspace_name}"
    }
    output {
        String upsert_entities_response = stdout()
    }
    runtime {
        docker: container
        cpu: select_first([runtime_set_cpu, 1])
        gpu: false
        memory: select_first([runtime_set_memory, 2]) + " GB"
        disks: "local-disk " + select_first([runtime_set_disk, 10]) + if use_ssd then " SSD" else " HDD"
        maxRetries: select_first([runtime_set_max_retries, 0])
        preemptible: select_first([runtime_set_preemptible_tries, 5])
        returnCodes: 0
     }
}
