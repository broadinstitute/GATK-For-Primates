version development

## Define the primary workflow
workflow GATKForAnimalsGermlineSNPsIndels_GATK4 {

    ####################################################################
    ## Initial workflow setup and configuration
    ####################################################################

    ## Workflow inputs
    input {

        ## Collect workflow options from input JSON
        String mode # Options: initial repeat final
        Boolean bwamem2 = false # Options: true / false to use either bwa (as bwa mem) or bwamem2 (as bwamem2 mem)
        Boolean flowcell_patterned = false # Options: true / false to use either bwa (as bwa mem) or bwamem2 (as bwamem2 mem)
        Int? merge_contigs_into_num_partitions # Option to speed up GenomicsDBImport with large number of contigs
        File? truth_set_SNPs # Option to provide a SNP truth set for VQSR (training set is produced via hard filtering)
        File? truth_set_INDELs # Option to provide an INDEL truth set for VQSR (training set is produced via hard filtering)

        ## Collect reference input file
        File ref

        ## Define arrays from input JSON; structural definitions are at the end of workflow
        Array[sampleInfo] sampleList
        Array[scatterInfo] scatterList
        
        ## Placeholders for when I add runtime/docker info
        ## Is it just me or has the 'Genomes in the Cloud' docker from Broad not been updated
        ## in three years? https://hub.docker.com/r/broadinstitute/genomes-in-the-cloud
        ## So using this one instead for bwa and samtools
        String docker_image_gatk = "broadinstitute/gatk:4.2.0.0"
        String docker_image_bwa_and_samtools = "bioslimcontainers/bwa-samtools:bwa-0.7.17_samtools-1.11"
        String docker_image_python = "python:latest"

    }

    ## Determine correct reference index files depending on version of bwa is being used
    File ref_dict = sub(ref, ".fa", ".dict")
    Array[File] ref_idxs = if bwamem2 then prefix(ref + ".", ["fai", "amb", "ann", "pac", "0123","bwt.2bit.64"]) else prefix(ref + ".", ["fai", "amb", "ann", "pac", "bwt", "sa"])

    ## Set correct command for executing bwa depending on which version is being used
    String execute_bwa = if bwamem2 then "bwa-mem2" else "bwa"

    ## Set configuration Booleans as true or false
    Boolean input_is_fastq = mode == "initial"
    Boolean truth_SNPs = defined(truth_set_SNPs)
    Boolean truth_INDELs = defined(truth_set_INDELs)

    ####################################################################
    ## Initial scatters to collect variables needed later scattering
    ####################################################################

    ## Scatter over scatters to get list of scatter names and intervals
    ## All this is necessary to sub-scatter later on
    scatter (sc in scatterList) {
        call collectScatters {
            input:
                scatterName = sc.name,
                scatterIntervals = sc.intervals,
        }
    }

    ## Define arrays comprising all scatter names and intervals
    ## Array[String] allScatterNames = collectScatters.allScatterNames
    ## Array[String] allScatterIntervals = collectScatters.allScatterIntervals

    ## Scatter over scatters to get list of group and sample names
    ## All this is necessary to sub-scatter later on and for counting number of samples
    scatter (sample in sampleList) {
        call collectAllGroupNames {
        input:
            groupName = sample.group,
            sampleName = sample.name,
        }
    }
    ## Array[String] allGroupsForScattering = collectAllGroupNames.allGroupNames
    Array[String] allSampleNames = collectAllGroupNames.allSampleNames


    ####################################################################
    ## MODE = INITIAL: Map FASTQ files to BAM
    ####################################################################
    
    if (input_is_fastq && mode == "initial") {

        ## Scatter over samples
        scatter (sample in sampleList) {

            ## Map paired-end reads into BAM file
            call FastqToBwaMem {
            input:
                ref = ref,
                ref_dict = ref_dict,
                ref_idxs = ref_idxs,
                sampleName = sample.name,
                rgID = sample.rgID,
                rgSM = sample.rgSM,
                rgLB = sample.rgLB,
                R1 = sample.R1,
                R2 = sample.R2,
                execute_bwa = execute_bwa,
                docker_image = docker_image_bwa_and_samtools,
            }

            ## Run MarkDuplicatesSpark
            ## Outputs BAM file sorted by co-ordinate
            call MarkDuplicatesSpark {
            input:
                sampleName = sample.name,
                input_bam = FastqToBwaMem.output_bam,
                flowcell_patterned = flowcell_patterned,
                docker_image = docker_image_gatk,
            }

            ## Run SetNmMdAndUqTags
            ## Outputs BAM file with fixed tags
            call SortAndFixTags {
            input:
                ref = ref,
                ref_dict = ref_dict,
                ref_idxs = ref_idxs,
                sampleName = sample.name,
                sampleGroup = sample.group,                
                input_bam = MarkDuplicatesSpark.output_bam,
                docker_image = docker_image_gatk,
            }

        }

        ## Generate a TSV (not really a 'map' per se) of the newly created BAM files
        call GenerateMapOfNewBAMs {
        input:
            new_bams = SortAndFixTags.output_bam,
            new_bam_groups = SortAndFixTags.output_group,
            new_bam_names = SortAndFixTags.output_name,
            new_bam_indexes = SortAndFixTags.output_bam_index,
        }

        ## Get array of all group names for scattering over later
        #Array[String] allGroupsForScattering = GenerateMapOfNewBAMs.allGroupsForScattering

        ## Coerce the file from File? into File because Cromwell cannot accommodate
        ## using files that may not exist in ternary operators (i.e. if-then-else)
        call coerceFile as coerceMapOfNewBAMs {
        input:
            input_file = GenerateMapOfNewBAMs.map,
        }

        ## Read the new TSV into a map
        Array[Array[String]] map_of_new_bams = read_tsv(coerceMapOfNewBAMs.coerced)

    }

    ####################################################################
    ## MODE = REPEAT/FINAL: Use existing BAM files
    ####################################################################

    if (!input_is_fastq) {

        ## Collect info on each sample from the input JSON file
        scatter (sample in sampleList) {
            call collectExistingBamInfo {
                input:
                    sampleName = sample.name,
                    sampleGroup = sample.group,
                    sampleBAM = sample.bam,
                    sampleBAI = sample.bam_index,
            }
        }

        ## Pull all this information together into arrays
        Array[String] existingSamples = collectExistingBamInfo.allNames
        Array[String] existingGroups = collectExistingBamInfo.allGroups
        Array[String] existingBAMs = collectExistingBamInfo.allBams
        Array[String] existingBAIs = collectExistingBamInfo.allBais

        ## Use the arrays to generate a TSV (not really a 'map' per se) of the existing BAM files
        call GenerateMapOfExistingBAMs {
        input:
            existing_bam_names = existingSamples,
            existing_bam_groups = existingGroups,
            existing_bams = existingBAMs,
            existing_bam_indexes = existingBAIs,
        }

        ## Coerce the file from File? into File because Cromwell cannot accommodate
        ## using files that may not exist in ternary operators (i.e. if-then-else)
        call coerceFile as coerceMapOfExistingBAMs {
        input:
            input_file = GenerateMapOfExistingBAMs.map,
        }

        ## Read the new TSV into a map
        Array[Array[String]] map_of_existing_bams = read_tsv(coerceMapOfExistingBAMs.coerced)

    }

    ####################################################################
    ## MODE = ANY: Work out which BAM files to use hereon out
    ####################################################################

    ## Work out which map of BAM files to use; either new or existing
    ## This has to be passed to a task because Cromwell cannot accommodate
    ## using files that may not exist in ternary operators (i.e. if-then-else)
    call FinalizeInputs {
    input:
        map_of_new_bams = map_of_new_bams,
        map_of_existing_bams = map_of_existing_bams,
    }

    ####################################################################
    ## MODE = ANY: Run HaplotypeCaller per sample per scatter interval
    ####################################################################

    ## Now we scatter over the user-defined scatter intervals
    scatter (sc in scatterList) {

        ## Then we scatter over each input BAM file
        scatter (col in FinalizeInputs.finalized_inputs) {

            ## Call HaplotypeCaller
            ## This calls haplotypes for the current sample at the current scatter
            call HaplotypeCaller {
            input:
                ref = ref,
                ref_dict = ref_dict,
                ref_idxs = ref_idxs,
                input_bam = col[2],
                input_bam_index = col[3],
                gvcf_basename = col[0] + ".haplotypes." + col[1] + "." + sc.name,
                sampleName = col[0],
                sampleGroup = col[1],
                scatterName = sc.name,
                scatterIntervals = sub(sc.intervals, " ", " -L "),
                docker_image = docker_image_gatk,
            }

        }
        ## We've just stopped scattering over samples but are still scattering
        ## over the user-defined scatter interval.

        ## Now we scatter over each group as pair.left and send an array of all member samples as pair.right
        scatter (pair in as_pairs(collect_by_key(zip(HaplotypeCaller.output_gvcf_group,HaplotypeCaller.output_gvcf)))) {

            ## Now we import all gVCFS per group into a GenomicsDB
            ## There is a bit of a bottleneck here; because we used pair and collect_by_key we cannot send the
            ## .tbi file with the .gvcf.gz -- so this is created again in the task.
            call GenomicsDBImport {
            input:
                database_name = "initial_variant_calls_" + pair.left + "_" + sc.name,
                input_gvcfs = pair.right,
                scatterIntervals = sc.intervals,
                docker_image = docker_image_gatk,
                merge_contigs_into_num_partitions = merge_contigs_into_num_partitions,
            }

            call GenotypeGenomicsDB {
            input:
                ref = ref,
                ref_dict = ref_dict,
                ref_idxs = ref_idxs,
                input_genomicsdb = GenomicsDBImport.output_genomicsdb,
                database_name = "initial_variant_calls_" + pair.left + "_" + sc.name,
                scatterIntervals = sc.intervals,
                docker_image = docker_image_gatk,
                groupName = pair.left,
            }

            ## This only really needs to be done if it's the final run
            ## or when there is no defined SNP truth set, but WDL
            ## makes this difficult because it doesn't like conditionals
            ## -- and it doesn't take very long anyway
            call HardFilterSNPs {
            input:
                ref = ref,
                ref_dict = ref_dict,
                ref_idxs = ref_idxs,
                input_genotypes = GenotypeGenomicsDB.output_genotypes,
                input_genotypes_index = GenotypeGenomicsDB.output_genotypes_index,
                groupName = pair.left,
                scatterName = sc.name,
                scatterIntervals = sc.intervals,
                docker_image = docker_image_gatk,
            }

            call HardFilterINDELs {
            input:
                ref = ref,
                ref_dict = ref_dict,
                ref_idxs = ref_idxs,
                input_genotypes = GenotypeGenomicsDB.output_genotypes,
                input_genotypes_index = GenotypeGenomicsDB.output_genotypes_index,
                groupName = pair.left,
                scatterName = sc.name,
                scatterIntervals = sc.intervals,
                docker_image = docker_image_gatk,
            }

        } ## Stop scattering over group/sample

    } ## Stop scattering over scatter list

    ####################################################################
    ## MODE = INITIAL/REPEAT: Base recalibration
    ####################################################################

    if (mode == "initial" || mode == "repeat") {

        ## Gather hard filtered SNPs
        scatter (pair in as_pairs(collect_by_key(zip(flatten(HardFilterSNPs.output_groupName),flatten(HardFilterSNPs.output_filtered_SNPs))))) {

            call GatherHardFilteredVcfs as GatherSNPs {
            input:
                input_filtered = pair.right,
                ref_dict = ref_dict,
                docker_image = docker_image_gatk,
                groupName = pair.left,
                type = "SNP",
            }
        }

        scatter (pair in as_pairs(collect_by_key(zip(flatten(HardFilterINDELs.output_groupName),flatten(HardFilterINDELs.output_filtered_INDELs))))) {

            call GatherHardFilteredVcfs as GatherINDELs {
            input:
                groupName = pair.left,
                input_filtered = pair.right,
                ref_dict = ref_dict,
                docker_image = docker_image_gatk,
                type = "INDEL",
            }
        }

        ## Then we scatter over each input BAM file to recalibrate
        scatter (col in FinalizeInputs.finalized_inputs) {

            ## Call BaseRecalibrator etc. for each BAM file
            call BaseRecalibrator {
            input:
                    ref = ref,
                    ref_dict = ref_dict,
                    ref_idxs = ref_idxs,
                    input_bam = col[2],
                    input_bam_index = col[3],
                    input_SNP_sites = GatherSNPs.output_filtered_sites_only,
                    input_SNP_sites_indexes = GatherSNPs.output_filtered_sites_only_index,
                    input_INDEL_sites = GatherINDELs.output_filtered_sites_only,
                    input_INDEL_sites_indexes = GatherINDELs.output_filtered_sites_only_index,
                    sampleName = col[0],
                    docker_image = docker_image_gatk,
            }
        }

    }

    ## Workflow is now complete if mode is initial or repeat
    
    ####################################################################
    ## MODE = FINAL: Variant recalibration/filtration
    ####################################################################

    if (mode == "final" && truth_SNPs) {

        scatter (pair in as_pairs(collect_by_key(zip(flatten(GenotypeGenomicsDB.output_groupName),flatten(GenotypeGenomicsDB.output_genotypes))))) {

            call GatherUnfilteredVcfs as GatherUnfilteredSNPs {
            input:
                groupName = pair.left,
                input_unfiltered = pair.right,
                type = "SNP",
                ref = ref,
                ref_dict = ref_dict,
                ref_idxs = ref_idxs,
                docker_image = docker_image_gatk,
            }

        }

        call VariantRecalibrator as RecalibrateSNPs {
        input:
            ref = ref,
            ref_dict = ref_dict,
            ref_idxs = ref_idxs,
            truth_set = truth_set_SNPs,
            input_sites = GatherUnfilteredSNPs.output_unfiltered_sites_only,
            input_sites_indexes = GatherUnfilteredSNPs.output_unfiltered_sites_only_index,
            type = "SNP",
            docker_image = docker_image_gatk,
        }

    } 

    if (mode == "final" && truth_INDELs) {

        scatter (pair in as_pairs(collect_by_key(zip(flatten(GenotypeGenomicsDB.output_groupName),flatten(GenotypeGenomicsDB.output_genotypes))))) {

            call GatherUnfilteredVcfs as GatherUnfilteredINDELs {
            input:
                groupName = pair.left,
                input_unfiltered = pair.right,
                type = "INDEL",
                ref = ref,
                ref_dict = ref_dict,
                ref_idxs = ref_idxs,
                docker_image = docker_image_gatk,

            }
        } 

        call VariantRecalibrator as RecalibrateINDELs {
        input:
            ref = ref,
            ref_dict = ref_dict,
            ref_idxs = ref_idxs,
            truth_set = truth_set_INDELs,
            input_sites = GatherUnfilteredINDELs.output_unfiltered_sites_only,
            input_sites_indexes = GatherUnfilteredINDELs.output_unfiltered_sites_only_index,
            type = "INDEL",
            docker_image = docker_image_gatk,
        }

    }



####################################################################
   ## MODE = FINAL: Perform final round of variant calls
####################################################################

    call ProduceFinalCallset {
    input:
        ref_dict = ref_dict,
        SNPs_hard_filtered = GatherSNPs.output_filtered_sites_only,
        SNPs_hard_filtered_index = GatherSNPs.output_filtered_sites_only_index,
        INDELs_hard_filtered = GatherINDELs.output_filtered_sites_only,
        INDELs_hard_filtered_index = GatherINDELs.output_filtered_sites_only_index,
        SNPs_recalibrated = RecalibrateSNPs.output_recalibrated_sites_only,
        SNPs_recalibrated_index = RecalibrateSNPs.output_recalibrated_sites_only_index,
        INDELs_recalibrated = RecalibrateINDELs.output_recalibrated_sites_only,
        INDELs_recalibrated_index = RecalibrateINDELs.output_recalibrated_sites_only_index,
    }


    scatter (col in FinalizeInputs.finalized_inputs) {

        ## Call HaplotypeCaller
        ## This calls haplotypes within 100 bp -ip of all the final callset sites
        call HaplotypeCaller as HaplotypeCallerFinal {
        input:
            ref = ref,
            ref_dict = ref_dict,
            ref_idxs = ref_idxs,
            input_bam = col[2],
            input_bam_index = col[3],
            gvcf_basename = col[0] + ".haplotypes",
            sampleName = col[0],
            FinalCallset = ProduceFinalCallset.FinalCallset,
            FinalCallset_index = ProduceFinalCallset.FinalCallset_index,
            docker_image = docker_image_gatk,
        }

    }

    call GenomicsDBImportFinal {
    input:
        input_gvcfs = HaplotypeCallerFinal.output_gvcf,
        input_gvcfs_indexes = HaplotypeCallerFinal.output_gvcf_index,
        database_name = "final_variant_calls",
        FinalCallset = ProduceFinalCallset.FinalCallset,
        FinalCallset_index = ProduceFinalCallset.FinalCallset_index,
        docker_image = docker_image_gatk,
    }

    call GenotypeGenomicsDB as GenotypeGenomicsDBFinal {
    input:
        ref = ref,
        ref_dict = ref_dict,
        ref_idxs = ref_idxs,
        input_genomicsdb = GenomicsDBImportFinal.output_genomicsdb,
        database_name = "final_variant_calls",
        FinalCallset = ProduceFinalCallset.FinalCallset,
        FinalCallset_index = ProduceFinalCallset.FinalCallset_index,        
        docker_image = docker_image_gatk,
    }

} ## End workflow

####################################################################
## Define tasks in the order that they are called
####################################################################

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

task FastqToBwaMem {
    input {
        File ref
        File ref_dict
        Array[File]+ ref_idxs
        String sampleName
        String? rgID
        String? rgSM
        String? rgLB
        File? R1
        File? R2
        String execute_bwa
        String docker_image
    }
    command {
        ${execute_bwa} mem -K 100000000 -Y -v 3 -R "@RG\tID:${rgID}\tPL:ILLUMINA\tLB:${rgLB}\tSM:${rgSM}" ${ref} ${R1} ${R2} | samtools sort -n -o ~{sampleName}.bam -

        
    }
    output {
        File output_bam = "~{sampleName}.bam"
        
    }
    runtime {
        docker: docker_image
    }
}

task MarkDuplicatesSpark {
    input {
        String docker_image
        File input_bam
        String sampleName
        Boolean flowcell_patterned
    }
    String pixel_distance = if flowcell_patterned then "2500" else "100"
        command {
        gatk \
        MarkDuplicatesSpark \
        -I ~{input_bam} \
        --optical-duplicate-pixel-distance ~{pixel_distance} \
        -O ~{sampleName}.dedup.bam \
        -M ~{sampleName}.metrics.txt
    }
    output {
        File output_bam = "~{sampleName}.dedup.bam"
        File duplication_metrics = "~{sampleName}.metrics.txt"
    }
    #runtime {
    #    docker: docker_image
    #}
}

task SortAndFixTags {
    input {
        File ref
        File ref_dict
        Array[File]+ ref_idxs
        String docker_image
        File input_bam
        String sampleName
        String sampleGroup
    }
    command {
        gatk \
        SetNmMdAndUqTags \
        -I ~{input_bam} \
        -R ~{ref} \
        -O ~{sampleName}.dedup.tagged.bam \
        --CREATE_INDEX true
    }
    output {
        File output_bam = "~{sampleName}.dedup.tagged.bam"
        File output_bam_index = "~{sampleName}.dedup.tagged.bai"
        String output_name = "~{sampleName}"        
        String output_group = "~{sampleGroup}"
    }
    runtime {
        docker: docker_image
    }
}

task GenerateMapOfNewBAMs {
## Adapted from the generate-sample-map.wdl (C) Broad Institute 2020
## https://github.com/gatk-workflows/utility-wdls/blob/main/generate-sample-map.wdl
    input {
        Array[String] new_bams
        Array[String] new_bam_indexes
        Array[String] new_bam_names
        Array[String] new_bam_groups
        String docker = "python:latest"
    }
    command <<<
    set -oe pipefail
    
    python << CODE
    new_bams = ['~{sep="','" new_bams}']
    new_bam_indexes = ['~{sep="','" new_bam_indexes}']
    new_bam_names = ['~{sep="','" new_bam_names}']
    new_bam_groups = ['~{sep="','" new_bam_groups}']

    if len(new_bams)!= len(new_bam_indexes) != len(new_bam_names) != len(new_bam_groups):
      print("Numbers of input variables are not equal")
      exit(1)

    with open("new_bams_map_file.txt", "w") as fi:
      for i in range(len(new_bams)):
        fi.write(new_bam_names[i] + "\t" + new_bam_groups[i] + "\t" + new_bams[i] + "\t" + new_bam_indexes[i] + "\n")
    CODE
    >>>

    runtime {
        docker: docker
    }

    output {
        File map = "new_bams_map_file.txt"
        #Array[String] bam_names = "~{new_bams}"
        #Array[String] bam_groups = "~{new_bam_groups}"
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

task collectExistingBamInfo {
    input {
        String sampleName
        String sampleGroup
        File? sampleBAM
        File? sampleBAI
    }
    command <<<
        ## The only purpose of this it coerce the ? inputs to outputs
    >>>
    output {
        String allNames = sampleName
        String allGroups = sampleGroup
        String allBams = sampleBAM
        String allBais = sampleBAI
    }
}


task GenerateMapOfExistingBAMs {
## Adapted from the generate-sample-map.wdl (C) Broad Institute 2020
## https://github.com/gatk-workflows/utility-wdls/blob/main/generate-sample-map.wdl
    input {
        Array[String] existing_bams
        Array[String] existing_bam_indexes
        Array[String] existing_bam_names
        Array[String] existing_bam_groups
        String docker = "python:latest"
    }
    command <<<
    set -oe pipefail
    
    python << CODE
    existing_bams = ['~{sep="','" existing_bams}']
    existing_bam_indexes = ['~{sep="','" existing_bam_indexes}']
    existing_bam_names = ['~{sep="','" existing_bam_names}']
    existing_bam_groups = ['~{sep="','" existing_bam_groups}']

    if len(existing_bams)!= len(existing_bam_indexes) != len(existing_bam_names) != len(existing_bam_groups):
      print("Numbers of input variables are not equal")
      exit(1)

    with open("existing_bams_map_file.txt", "w") as fi:
      for i in range(len(existing_bams)):
        fi.write(existing_bam_names[i] + "\t" + existing_bam_groups[i] + "\t" + existing_bams[i] + "\t" + existing_bam_indexes[i] + "\n")
    CODE
    >>>

    runtime {
        docker: docker
    }

    output {
        File map = "existing_bams_map_file.txt"
        #Array[String] bam_names = "~{existing_bams}"
        #Array[String] bam_groups = "~{existing_bam_groups}"
    }

}

task FinalizeInputs {
    input {
        Array[Array[String]]? map_of_new_bams
        Array[Array[String]]? map_of_existing_bams
    }
    Array[Array[String]]? map_to_use = if defined(map_of_new_bams) then map_of_new_bams else map_of_existing_bams
    command <<<
        ## The only purpose of this task is to create the Array[Array[String]] map_to_use
    >>>
    output {
        Array[Array[String]] finalized_inputs = map_to_use
    }
}

task HaplotypeCaller {
    input {
        File ref
        File ref_dict
        Array[File]+ ref_idxs
        String docker_image
        File? input_bam
        File? input_bam_index
        String gvcf_basename
        String sampleName
        String? sampleGroup
        String? scatterName
        String? scatterIntervals
        File? FinalCallset
        File? FinalCallset_index
    }
    String padding = if defined(FinalCallset) then "--interval-padding 100 \\" else " "

    command {
        gatk \
        HaplotypeCaller \
        -I ~{input_bam} \
        -R ~{ref} \
        -O ~{gvcf_basename}.g.vcf.gz \
        -L ~{scatterIntervals} ~{FinalCallset} \
        -ERC GVCF \
        -G StandardAnnotation \
        -G AS_StandardAnnotation \
        -G StandardHCAnnotation \
        -GQB 10 -GQB 20 -GQB 30 -GQB 40 -GQB 50 -GQB 60 -GQB 70 -GQB 80 -GQB 90 \
        ~{padding}
    }
    output {
        File output_gvcf = "~{gvcf_basename}" + ".g.vcf.gz"
        File output_gvcf_index = "~{gvcf_basename}" + ".g.vcf.gz.tbi"
        String output_gvcf_name = "~{sampleName}"
        String output_gvcf_scattername = "~{scatterName}"
        String output_gvcf_group = "~{sampleGroup}"
    }
    #runtime {
        #docker: docker_image
    #}
}

task GenomicsDBImport {
    input {
        Array[File] input_gvcfs
        String database_name
        String scatterIntervals
        String docker_image
        Int? merge_contigs_into_num_partitions
    }
    #Int merge_contigs_value = if defined(merge_contigs_into_num_partitions) then merge_contigs_into_num_partitions else "0"
    command <<<
    set -euo pipefail

        gatk IndexFeatureFile -I ~{sep=' && gatk IndexFeatureFile -I ' input_gvcfs}

        gatk \
        GenomicsDBImport \
        -V ~{sep=' -V ' input_gvcfs} \
        --genomicsdb-workspace-path ~{database_name} \
        -L ~{scatterIntervals}

        tar -cf ~{database_name}.tar ~{database_name}
    >>>
    output {
        File output_genomicsdb = "~{database_name}.tar"
    }
    #runtime {
    #    docker: docker_image
    #}
}

task GenotypeGenomicsDB {
    input {
        File ref
        File ref_dict
        Array[File]+ ref_idxs
        String database_name
        File input_genomicsdb
        File? FinalCallset
        File? FinalCallset_index
        String? scatterIntervals
        String docker_image
        String? groupName
    }
    String final_options = if defined(FinalCallset) then "--interval-padding 100 --include-non-variant-sites\\" else " "
    command <<<
    set -euo pipefail

        tar -xf ~{input_genomicsdb}

        gatk \
        GenotypeGVCFs \
        -R ~{ref} \
        -V gendb://~{database_name} \
        -L ~{scatterIntervals} ~{FinalCallset} \
        -O ~{database_name}.vcf.gz \
        -G StandardAnnotation -G AS_StandardAnnotation \
        ~{final_options}
        
    >>>
    output {
        String output_groupName = "~{groupName}"
        File output_genotypes = "~{database_name}.vcf.gz"
        File output_genotypes_index = "~{database_name}.vcf.gz.tbi"
    }
}

## Need to add AS annotation flags here
task HardFilterSNPs {
    input {
        File ref
        File ref_dict
        Array[File]+ ref_idxs
        File input_genotypes
        File input_genotypes_index
        String groupName
        String scatterName
        String scatterIntervals
        String docker_image
    }
    command {

        gatk \
        SelectVariants \
        -R ~{ref} \
        -V ~{input_genotypes} \
        --select-type-to-include SNP \
        -O ~{groupName}_~{scatterName}_SNPs_unfiltered.vcf.gz

        gatk \
        VariantFiltration \
        -R ~{ref} \
        -V ~{groupName}_~{scatterName}_SNPs_unfiltered.vcf.gz \
        -O ~{groupName}_~{scatterName}_SNPs_unfiltered_and_filtered.vcf.gz \
        -filter "QD < 2.0" --filter-name "QD2" \
        -filter "QUAL < 30.0" --filter-name "QUAL30" \
        -filter "SOR > 3.0" --filter-name "SOR3" \
        -filter "FS > 60.0" --filter-name "FS60" \
        -filter "MQ < 40.0" --filter-name "MQ40" \
        -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
        -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
        --verbosity ERROR

        gatk \
        SelectVariants \
        -R ~{ref} \
        -V ~{groupName}_~{scatterName}_SNPs_unfiltered_and_filtered.vcf.gz \
        -O ~{groupName}_~{scatterName}_SNPs_filtered.vcf.gz \
        --exclude-filtered
        
    }
    output {
        File output_filtered_SNPs = "~{groupName}_~{scatterName}_SNPs_filtered.vcf.gz"
        File output_filtered_SNPs_index = "~{groupName}_~{scatterName}_SNPs_filtered.vcf.gz.tbi"
        String output_groupName = "~{groupName}"
    }
}

## Need to add AS annotation flags here
task HardFilterINDELs {
    input {
        File ref
        File ref_dict
        Array[File]+ ref_idxs
        File input_genotypes
        File input_genotypes_index
        String groupName
        String scatterName
        String scatterIntervals
        String docker_image
    }
    command {

        gatk \
        SelectVariants \
        -R ~{ref} \
        -V ~{input_genotypes} \
        --select-type-to-include MIXED \
        -O ~{groupName}_~{scatterName}_INDELs_unfiltered.vcf.gz

        gatk \
        VariantFiltration \
        -R ~{ref} \
        -V ~{groupName}_~{scatterName}_INDELs_unfiltered.vcf.gz \
        -O ~{groupName}_~{scatterName}_INDELs_unfiltered_and_filtered.vcf.gz \
        -filter "QD < 2.0" --filter-name "QD2" \
        -filter "QUAL < 30.0" --filter-name "QUAL30" \
        -filter "FS > 200.0" --filter-name "FS200" \
        -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
        --verbosity ERROR

        gatk \
        SelectVariants \
        -R ~{ref} \
        -V ~{groupName}_~{scatterName}_INDELs_unfiltered_and_filtered.vcf.gz \
        -O ~{groupName}_~{scatterName}_INDELs_filtered.vcf.gz \
        --exclude-filtered
        
    }
    output {
        File output_filtered_INDELs = "~{groupName}_~{scatterName}_INDELs_filtered.vcf.gz"
        File output_filtered_INDELs_index = "~{groupName}_~{scatterName}_INDELs_filtered.vcf.gz.tbi"
        String output_groupName = "~{groupName}"
    }
}

task GatherHardFilteredVcfs {
    input {
        Array[File] input_filtered
        File ref_dict
        String groupName
        String type
        String docker_image
    }
    command {
    
        ## You might not need this is GatherVcfs doesn't require indexes
        ## Commented out to try it
        ## gatk IndexFeatureFile -I ~{sep=' && gatk IndexFeatureFile -I ' input_filtered}

            gatk \
            GatherVcfs \
            -I ~{sep=" -I " input_filtered} \
            -O ~{groupName}_~{type}_filtered.vcf.gz \
            --REORDER_INPUT_BY_FIRST_VARIANT

            gatk MakeSitesOnlyVcf \
            -I ~{groupName}_~{type}_filtered.vcf.gz \
            -O ~{groupName}_~{type}_filtered_sites_only_to_sort.vcf.gz

            ## This may not be necessary assuming the -RI flag in GatherVcfs achieves this
            ## Check this and revise if necessary
            gatk SortVcf \
            -I ~{groupName}_~{type}_filtered_sites_only_to_sort.vcf.gz \
            -O ~{groupName}_~{type}_filtered_sites_only.vcf.gz \
            -SD ~{ref_dict}

    }
    output {
        File output_filtered_sites_only = "~{groupName}_~{type}_filtered_sites_only.vcf.gz"
        File output_filtered_sites_only_index = "~{groupName}_~{type}_filtered_sites_only.vcf.gz.tbi"
    }
}

task BaseRecalibrator {
    input {
        File ref
        File ref_dict
        Array[File]+ ref_idxs
        File input_bam
        File input_bam_index
        String sampleName
        Array[File] input_SNP_sites
        Array[File] input_SNP_sites_indexes
        Array[File] input_INDEL_sites
        Array[File] input_INDEL_sites_indexes
        String docker_image
    }
    command <<<
        set -euo pipefail
        
        #GENERATE BEFORE TABLE FROM THE UNRECALIBRATED BAM FILE USING THE HARD FILTERED CALLS FROM THE CURRENT/LATEST PASS
        gatk \
        BaseRecalibrator \
        -I ~{input_bam} \
        -R ~{ref} \
        --known-sites ~{sep=" --known-sites " input_SNP_sites} \
        --known-sites ~{sep=" --known-sites " input_INDEL_sites} \
        -O ~{sampleName}.before.table

        #APPLY BQSR TO ADJUST THE QUALITY SCORES
        gatk \
        ApplyBQSR \
        -I ~{input_bam} \
        -R ~{ref} \
        --bqsr-recal-file ~{sampleName}.before.table \
        -O ~{sampleName}_recalibrated.bam

        #GENERATE AFTER TABLE
        gatk \
        BaseRecalibrator \
        -I ~{sampleName}_recalibrated.bam \
        -R ~{ref} \
        --known-sites ~{sep=" --known-sites " input_SNP_sites} \
        --known-sites ~{sep=" --known-sites " input_INDEL_sites} \
        -O ~{sampleName}.after.table

        #GENERATE PLOTS
        #gatk \
        #AnalyzeCovariates \
        #-before ~{sampleName}.before.table \
        #-after ~{sampleName}.after.table \
        #-plots ~{sampleName}.AnalyzeCovariates_plots.pdf \
        #--use-jdk-deflater \
        #--use-jdk-inflater
    >>>
    output {
        File table_before = "~{sampleName}.before.table"
        File table_after = "~{sampleName}.after.table"
        #File plots = "~{sampleName}.AnalyzeCovariates_plots.pdf"
        File recalibrated_bam = "~{sampleName}_recalibrated.bam"
        File recalibrated_bam_index = "~{sampleName}_recalibrated.bai"
    }
}

task GatherUnfilteredVcfs {
    input {
        Array[File] input_unfiltered
        File ref
        File ref_dict
        Array[File]+ ref_idxs
        String groupName
        String type
        String docker_image
    }
    String type_corrected = if type == "INDEL" then "MIXED" else "SNP"
    command {
    
            gatk \
            GatherVcfs \
            -I ~{sep="-I " input_unfiltered} \
            -O ~{groupName}_unfiltered.vcf.gz \
            --REORDER_INPUT_BY_FIRST_VARIANT

            gatk MakeSitesOnlyVcf \
            -I ~{groupName}_unfiltered.vcf.gz \
            -O ~{groupName}_unfiltered_sites_only_to_sort.vcf.gz

            ## This may not be necessary assuming the -RI flag in GatherVcfs achieves this
            ## Check this and revise if necessary
            gatk SortVcf \
            -I ~{groupName}_unfiltered_sites_only_to_sort.vcf.gz \
            -O ~{groupName}_unfiltered_sites_only.vcf.gz \
            -SD ~{ref_dict}

            gatk \
            SelectVariants \
            -R ~{ref} \
            -V ~{groupName}_unfiltered_sites_only.vcf.gz \
            --select-type-to-include ~{type} \
            -O ~{groupName}_~{type}_unfiltered_sites_only.vcf.gz

    }
    output {
        File output_unfiltered_sites_only = "~{groupName}_~{type}_unfiltered_sites_only.vcf.gz"
        File output_unfiltered_sites_only_index = "~{groupName}_~{type}_unfiltered_sites_only.vcf.gz.tbi"
    }
}

task VariantRecalibrator {
    input {
        File ref
        File ref_dict
        Array[File]+ ref_idxs
        File? truth_set
        Array[File] input_sites
        Array[File] input_sites_indexes
        String type
        String docker_image
    }
    command <<<
        set -euo pipefail

        if (type == "SNP") {

            gatk \
            VariantRecalibrator \
            -V ~{sep=" -V " input_sites} \
            --trust-all-polymorphic \
            -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.8 -tranche 99.6 -tranche 99.5 -tranche 99.4 -tranche 99.3 -tranche 99.0 -tranche 98.0 -tranche 97.0 -tranche 90.0 \
            -an QD -an MQRankSum -an ReadPosRankSum -an FS -an MQ -an SOR -an DP \
            -mode SNP \
            --max-gaussians 6 \
            -resource:truth_set_snps,known=false,training=true,truth=true,prior=~{truth_set} \
            -resource:hard_filtered,known=false,training=true,truth=false,prior=~{input_sites} \
            -O cohort_snps.recal \
            --tranches-file cohort_snps.tranches \
            --use-allele-specific-annotations


            gatk \
            ApplyVQSR \
            -V ~{sep=" -V " input_sites} \
            --recal-file cohort_snps.recal \
            --tranches-file cohort_snps.tranches \
            --truth-sensitivity-filter-level 99.7 \
            --create-output-variant-index true \
            -mode SNP \
            --exclude-filtered \
            --use-allele-specific-annotations \
            -O Cohort_~{type}.recalibrated.vcf.gz

        }

        if (type == "INDEL") {

            gatk \
            VariantRecalibrator \
            -V ~{sep=" -V " input_sites} \
            --trust-all-polymorphic \
            -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 97.0 -tranche 96.0 -tranche 95.0 -tranche 94.0 -tranche 93.5 -tranche 93.0 -tranche 92.0 -tranche 91.0 -tranche 90.0 \
            -an FS -an ReadPosRankSum -an MQRankSum -an QD -an SOR -an DP \      
            -mode MIXED \
            --max-gaussians 4 \
            -resource:truth_set_indels,known=false,training=true,truth=true,prior=~{truth_set} \
            -resource:hard_filtered,known=false,training=true,truth=false,prior=~{input_sites} \
            -O cohort_indels.recal \
            --use-allele-specific-annotations \
            --tranches-file cohort_indels.tranches

            gatk \
            ApplyVQSR \
            -V ~{sep=" -V " input_sites} \
            --recal-file cohort_indels.recal \
            --tranches-file cohort_indels.tranches \
            --truth-sensitivity-filter-level 99.7 \
            --create-output-variant-index true \
            -mode MIXED \
            --exclude-filtered \
            --use-allele-specific-annotations \
            -O Cohort_~{type}.recalibrated.vcf.gz

        }

        gatk MakeSitesOnlyVcf \
        -I Cohort_~{type}.recalibrated.vcf.gz \
        -O Cohort_~{type}.recalibrated_sites_only.vcf.gz


    >>>
    output {
        File output_recalibrated_sites_only = "Cohort_~{type}.recalibrated_sites_only.vcf.gz"
        File output_recalibrated_sites_only_index = "Cohort_~{type}.recalibrated_sites_only.vcf.gz.tbi"
    }
}

task ProduceFinalCallset {
    input {
        Array[File]? SNPs_hard_filtered
        Array[File]? SNPs_hard_filtered_index
        Array[File]? INDELs_hard_filtered
        Array[File]? INDELs_hard_filtered_index
        File? SNPs_recalibrated
        File? SNPs_recalibrated_index
        File? INDELs_recalibrated
        File? INDELs_recalibrated_index
        File ref_dict
    }
    String a = if defined(SNPs_hard_filtered) then " -I " else ""
    String b = if defined(SNPs_recalibrated) then " -I " else ""
    String c = if defined(INDELs_hard_filtered) then " -I " else ""
    String d = if defined(INDELs_recalibrated) then " -I " else ""
    command <<<
        set -euo pipefail

        gatk MergeVcfs \
        ~{a} ~{sep=" -I " SNPs_hard_filtered} ~{b} ~{sep=" -I " SNPs_recalibrated} ~{c} ~{sep=" -I " INDELs_hard_filtered} ~{d} ~{sep=" -I " INDELs_recalibrated} \
        -O FinalCallset.vcf.gz

    >>>
    output {
        File FinalCallset = "FinalCallset.vcf.gz"
        File FinalCallset_index = "FinalCallset.vcf.gz.tbi"
    }
}


task GenomicsDBImportFinal {
    input {
        Array[File] input_gvcfs
        Array[File] input_gvcfs_indexes
        String database_name
        File FinalCallset
        File FinalCallset_index
        String docker_image
    }
    command <<<
    set -euo pipefail

        gatk \
        GenomicsDBImport \
        -V ~{sep=' -V ' input_gvcfs} \
        --genomicsdb-workspace-path ~{database_name} \
        -L ~{FinalCallset} \
        --interval-padding 100

        tar -cf ~{database_name}.tar ~{database_name}
    >>>
    output {
        File output_genomicsdb = "~{database_name}.tar"
    }
    #runtime {
    #    docker: docker_image
    #}
}

####################################################################
## Define necessary structs
####################################################################

struct sampleInfo {
    String name
    String? rgID
    String? rgLB
    String? rgSM
    String group
    File? R1
    File? R2
    File? bam
    File? bam_index
}

struct scatterInfo {
    String name
    String intervals
}
