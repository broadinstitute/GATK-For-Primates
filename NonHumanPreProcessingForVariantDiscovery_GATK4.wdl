version development

## Workflow definition
workflow NonHumanPreProcessingForVariantDiscovery_GATK4 {

    ## Workflow inputs
    input {

        ## Define arrays from input JSON; structural definitions are at end of workflow
        Array[sampleInfo] sampleList
        Array[scatterInfo] scatterList

        ## Define files from input JSON
        File ref_fasta
        File ref_fasta_index
        File ref_dict
        File ref_sa 
        File ref_ann
        File ref_bwt
        File ref_pac
        File ref_amb

        ## Define required strings from input JSON
        String merge_contigs_into_num_partitions
        
        ## Placeholders for when I add runtime/docker info
        String gatk_docker = "broadinstitute/gatk:4.2.0.0"
        ## Is it just me or has the 'Genomes in the Cloud' docker from Broad not been updated
        ## in three years? https://hub.docker.com/r/broadinstitute/genomes-in-the-cloud
        ## So using this one instead for bwa and samtools
        String bwa_and_samtools_docker = "bioslimcontainers/bwa-samtools:bwa-0.7.17_samtools-1.11"

    }

    ## Scatter over samples just to get group names
    scatter (sample in sampleList) {

        ## Define strings comprising all data for each sample record
        String sampleName = sample.name
        String sampleGroup = sample.group

        ## Get complete list of samples and groups
        ## This is needed for scattering over groups later on
        call collectGroupsAndSamples {
            input:
                sampleName = sampleName,
                sampleGroup = sampleGroup,
        }
    }

    ## Define arrays for group names and group members
    ## This facilitates scattering over groups
    Array[String] groupNames = collectGroupsAndSamples.allGroups
    Array[String] groupSamples = collectGroupsAndSamples.allNames

    ## Scatter over groups and group samples
    ## Note this uses the collect_by_key function only available in the development version of WDL
    ## Whoever decided to put this function into WDL is a genius who must be feted
    ## I spent two weeks trying to get this to work until I found the dev spec!
    scatter (pair in as_pairs(collect_by_key(zip(groupNames,groupSamples)))) {

        String currentGroup = pair.left

        ## Scatter over samples
        scatter (sample in sampleList) {

            ## Define strings comprising all data for each sample record
            String sampleName = sample.name
            String rgID = sample.rgID
            String rgLB = sample.rgLB
            String rgSM = sample.rgSM
            String sampleGroup = sample.group
            String R1 = sample.R1
            String R2 = sample.R2

            ## Perform these actions only for samples within the specified group
            ## I know this is a pretty backwards way of doing things but WDL is not
            ## very friendly for scattering over maps and arrays
            if (currentGroup == sampleGroup) {

                ## Run bwa mem to map FASTQs to reference with read groups
                ## Outputs BAM file sorted by query name, as preferred by MarkDuplicatesSpark
                call FastqToBwaMem {
                    input:
                        ref_fasta = ref_fasta,
                        ref_fasta_index = ref_fasta_index,
                        ref_dict = ref_dict,
                        ref_sa = ref_sa,
                        ref_ann = ref_ann,
                        ref_bwt = ref_bwt,
                        ref_pac = ref_pac,
                        ref_amb = ref_amb,
                        docker_image = bwa_and_samtools_docker,
                        sampleName = sampleName,
                        rgID = rgID,
                        rgSM = rgSM,
                        rgLB = rgLB,
                        R1 = R1,
                        R2 = R2,
                }

                ## Run MarkDuplicatesSpark
                ## Outputs BAM file sorted by co-ordinate
                call MarkDuplicatesSpark {
                    input:
                        input_bam = FastqToBwaMem.output_bam,
                        docker_image = gatk_docker,
                        sampleName = sampleName,
                }

                ## Run SetNmMdAndUqTags
                ## Outputs BAM file with fixed tags
                call SortAndFixTags {
                    input:
                       input_bam = MarkDuplicatesSpark.output_bam,
                        ref_dict = ref_dict,
                        ref_fasta = ref_fasta,
                        ref_fasta_index = ref_fasta_index,
                        docker_image = gatk_docker,
                        sampleName = sampleName,
                }

                ## Scatter by each interval
                scatter (sc in scatterList) {

                    String scatterName = sc.name
                    String scatterIntervals = sub(sc.intervals, " ", " -L ")
                    String gvcf_basename = sampleName + ".haplotypes." + sampleGroup + "." + scatterName + ".g.vcf.gz"

                    call HaplotypeCaller {
                        input:
                            input_bam = SortAndFixTags.output_bam,
                            input_bam_index = SortAndFixTags.output_bam_index,
                            gvcf_basename = gvcf_basename,
                            scatterIntervals = scatterIntervals,
                            ref_dict = ref_dict,
                            ref_fasta = ref_fasta,
                            ref_fasta_index = ref_fasta_index,
                            docker_image = gatk_docker,
                    }

                }


            }

        }


        ## Flatten horrific nested array into string
        ## I know this is the world's jankiest code
        ## I will re-visit it later
        ## Also, this is annoying because it takes every HaplotypeCaller
        ## output as input files to the new job, even though
        ## you really only need the ones that are for this particular group
        ## but it doesn't appear to be possible to use a variable
        ## inside a call/task input, e.g. HaplotypeCaller.out_gvcf is fine
        ## where out_gvcfs is defined as the task output, but
        ## HaplotypeCaller.{$var_output} for selected output is not
        ## possible in either task output or call input
        Array[Array[File]?] gVCFS = HaplotypeCaller.output_gvcf
        String gVCFS_string_ugly = "~{sep=' -V ' gVCFS}"
        String gVCFS_string_ugly2 = sub(gVCFS_string_ugly, "\\[", "")
        String gVCFS_string_ugly3 = sub(gVCFS_string_ugly2, "\\]", "")
        String gVCFS_string_ugly4 = sub(gVCFS_string_ugly3, "\"", "")
        String gVCFS_string_ugly5 = sub(gVCFS_string_ugly4, ", ", " -V ")

        scatter (sc in scatterList) {

            String scatterName = sc.name
            String scatterIntervals = sub(sc.intervals, " ", " -L ")

            call GenomicsDBImport {
            input:
                input_gvcfs = gVCFS,
                input_gvcfs_string = gVCFS_string_ugly5,
                groupName = currentGroup,
                scatterName = scatterName,
                scatterIntervals = scatterIntervals,
                docker_image = gatk_docker,
            }

            ## The bottleneck ends here because it only takes the required genomicsDB
            call GenotypeGenomicsDB {
            input:
                ref_fasta = ref_fasta,
                ref_fasta_index = ref_fasta_index,
                ref_dict = ref_dict,
                input_genomicsdb = GenomicsDBImport.output_genomicsdb,
                groupName = currentGroup,
                scatterName = scatterName,
                scatterIntervals = scatterIntervals,
                docker_image = gatk_docker,
            }




        }





    ## End scatter by group
    }

## End workflow
}


    



## Define tasks





task collectGroupsAndSamples {
    input {
        String sampleName
        String sampleGroup
    }
    command {

    }
    output {
        String allNames = sampleName
        String allGroups = sampleGroup
    }
}

task FastqToBwaMem {
    input {
        File ref_fasta
        File ref_fasta_index
        File ref_dict
        File ref_amb
        File ref_ann
        File ref_bwt
        File ref_pac
        File ref_sa
        String docker_image
        String sampleName
        String rgID
        String rgSM
        String rgLB
        File R1
        File R2
    }
    command {
        bwa mem -K 100000000 -Y -v 3 -R "@RG\tID:${rgID}\tPL:ILLUMINA\tLB:${rgLB}\tSM:${rgSM}" ${ref_fasta} ${R1} ${R2} | samtools sort -n -o ~{sampleName}.bam -
        
    }
    output {
        File output_bam = "~{sampleName}.bam"
        
    }
}

task MarkDuplicatesSpark {
    input {
        String docker_image
        File input_bam
        String sampleName
    }    
    command {
        gatk \
        MarkDuplicatesSpark \
        -I ~{input_bam} \
        -O ~{sampleName}.dedup.bam \
        -M ~{sampleName}.metrics.txt
    }
    output {
        File output_bam = "~{sampleName}.dedup.bam"
        File duplication_metrics = "~{sampleName}.metrics.txt"
  }
}

task SortAndFixTags {
    input {
        File ref_fasta
        File ref_fasta_index
        File ref_dict
        String docker_image
        File input_bam
        String sampleName
    }
    command {
        gatk \
        SetNmMdAndUqTags \
        -I ~{input_bam} \
        -R ~{ref_fasta} \
        -O ~{sampleName}.dedup.tagged.bam \
        --CREATE_INDEX true
    }
    output {
        File output_bam = "~{sampleName}.dedup.tagged.bam"
        File output_bam_index = "~{sampleName}.dedup.tagged.bai"
    }
}

task HaplotypeCaller {
    input {
        File ref_fasta
        File ref_fasta_index
        File ref_dict
        String docker_image
        File input_bam
        File input_bam_index
        String gvcf_basename
        String scatterIntervals
    }
    command {
        gatk \
        HaplotypeCaller \
        -I ~{input_bam} \
        -R ~{ref_fasta} \
        -O ~{gvcf_basename} \
        -L ~{scatterIntervals} \
        -ERC GVCF \
        -G StandardAnnotation \
        -G AS_StandardAnnotation \
        -G StandardHCAnnotation \
        -GQB 10 -GQB 20 -GQB 30 -GQB 40 -GQB 50 -GQB 60 -GQB 70 -GQB 80 -GQB 90
    }
    output {
        File output_gvcf = "~{gvcf_basename}"
    }
}


## Need to add the new merge flag here and optimize
task GenomicsDBImport {
    input {
        Array[Array[File]?] input_gvcfs
        String input_gvcfs_string
        String groupName
        String scatterName
        String scatterIntervals
        String docker_image
    }
    command <<<
    set -euo pipefail

        gatk \
        GenomicsDBImport \
        -V ~{sep=' -V ' input_gvcfs_string} \
        --genomicsdb-workspace-path initial_variant_calls_~{groupName}_~{scatterName} \
        -L ~{scatterIntervals}

        tar -cf initial_variant_calls_~{groupName}_~{scatterName}.tar initial_variant_calls_~{groupName}_~{scatterName}/
    >>>
    output {
        File output_genomicsdb = "initial_variant_calls_~{groupName}_~{scatterName}.tar"
    }
}

## Need to add AS annotation flags here
task GenotypeGenomicsDB {
    input {
        File ref_fasta
        File ref_fasta_index
        File ref_dict
        File input_genomicsdb
        String groupName
        String scatterName
        String scatterIntervals
        String docker_image
    }
    command <<<
    set -euo pipefail

        tar -xf initial_variant_calls_~{groupName}_~{scatterName}.tar

        gatk \
        GenotypeGVCFs \
        -R ~{ref_fasta} \
        -V gendb://initial_variant_calls_~{groupName}_~{scatterName} \
        -L ~{scatterIntervals} \
        -O initial_variant_calls_~{groupName}_~{scatterName}.vcf.gz
        
    >>>
    output {
        File output_genotypes = "initial_variant_calls_~{groupName}_~{scatterName}.vcf.gz"
    }
}


## Define necessary structs

struct sampleInfo {
    String name
    String rgID
    String rgLB
    String rgSM
    String group
    File R1
    File R2
}

struct scatterInfo {
    String name
    String intervals
}
