version development

## Copyright Broad Institute and Wisconsin National Primate Research Center,
## University of Wisconsin-Madison, 2021
## 
## Sub-workflow from the complete germline short variant discovery pipeline
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
import "../tasks/Alignment.wdl" as alignment
import "../tasks/ProcessBAMs.wdl" as processBAMs
import "../tasks/Utilities.wdl" as utilities

##########################################################################
## WORKFLOW DEFINITION
##########################################################################

workflow unmappedInputToAlignedBAM {

    input {
        File ref
        File ref_dict
        File ref_fai
        Array[File?]+ ref_idxs
        Array[sampleInfo]+ sampleList
        Boolean flowcell_patterned
        Boolean bwamem2
        String container_gatk
        String container_gitc
        String container_python
        String path_to_gitc
        String path_to_gitc_gatk
    }

    ##########################################################################
    ## Set bwa command line and parameters, depending on aligner
    ##########################################################################

    String bwa_commandline = "bwa mem -K 100000000 -p -v 3 -t 16 -Y"
    String bwamem2_commandline = "TBD"
    String execute_aligner = if bwamem2 then bwamem2_commandline else bwa_commandline

    ##########################################################################
    ## Get bwa version. This is really only needed if uBAMs are provided
    ##########################################################################

    if (!bwamem2) {
        call alignment.getBwaVersion as getBwaVersion {
            input:
                # Runtime options
                container = container_gitc,
                path_to_gitc = path_to_gitc,
        }
    }

    ##########################################################################
    ## Map either the paired FASTQ or uBAM files for each sample
    ##########################################################################

    ## Scatter over each sample
    scatter (sample in sampleList) {

        # If no uBAM is available, map the FASTQ files
        if (!defined(sample.unmapped_bam)) {

            call alignment.mapFromPairedFASTQ as mapFromPairedFASTQ{
                input:
                    ref = ref,
                    ref_dict = ref_dict,
                    ref_fai = ref_fai,
                    ref_idxs = ref_idxs,
                    sampleName = sample.name,
                    RG_ID = sample.RG_ID,
                    RG_SM = sample.RG_SM,
                    RG_LB = sample.RG_LB,
                    RG_PU = sample.RG_PU,
                    R1 = sample.R1,
                    R2 = sample.R2,
                    execute_aligner = execute_aligner,
                    # Runtime options
                    container = container_gitc,
                    path_to_gitc = path_to_gitc,
            }

        }

        # Otherwise, map the uBAM instead
        if (defined(sample.unmapped_bam)) {

            ## Stream unmapped BAM to FASTQ and map to BAM
            call alignment.mapFromUnmappedBAM as mapFromUnmappedBAM {
                input:
                    ref = ref,
                    ref_dict = ref_dict,
                    ref_fai = ref_fai,
                    ref_idxs = ref_idxs,
                    sampleName = sample.name,
                    unmapped_bam = sample.unmapped_bam,
                    execute_aligner = execute_aligner,
                    # Runtime
                    container = container_gitc,
                    path_to_gitc = path_to_gitc,
                    path_to_gitc_gatk = path_to_gitc_gatk,
            }

            ## Merge mapped and unmapped BAMs
            call alignment.mergeMappedAndUnmapped as mergeMappedAndUnmapped{
                input:
                    ref = ref,
                    ref_dict = ref_dict,
                    ref_fai = ref_fai,
                    sampleName = sample.name,
                    unmapped_bam = sample.unmapped_bam,
                    mapped_unmerged_bam = mapFromUnmappedBAM.output_mapped_unmerged_bam,
                    execute_aligner = execute_aligner,
                    bwa_version = getBwaVersion.bwa_version,
                    # Runtime
                    container = container_gatk,
            }

        }

        ## Run MarkDuplicatesSpark on output BAM; output deduped co-ordinate-sorted BAM
        call processBAMs.markDuplicatesSpark as markDuplicatesSpark {
            input:
                sampleName = sample.name,
                input_bam = select_first([mapFromPairedFASTQ.output_bam_mapped, mergeMappedAndUnmapped.output_bam_mapped]),
                flowcell_patterned = flowcell_patterned,
                # Runtime
                container = container_gatk,
        }

        ## Run SetNmMdAndUqTags, outputs BAM with fixed tags
        call processBAMs.sortAndFixTags as sortAndFixTags {
            input:
                ref = ref,
                ref_dict = ref_dict,
                ref_fai = ref_fai,
                sampleName = sample.name,
                sampleGroup = sample.taxon_group,                
                input_bam = markDuplicatesSpark.output_bam_dedup,
                # Runtime
                container = container_gatk,
        }
        
    }

    output {
        Array[String] sampleGroups = sortAndFixTags.output_sampleGroup
        Array[String] sampleNames = sortAndFixTags.output_sampleName
        Array[File] bams = sortAndFixTags.output_bam_dedup_tagged
        Array[File] bam_indexes = sortAndFixTags.output_bam_dedup_tagged_index
        Array[File] bam_dedup_metrics = markDuplicatesSpark.duplication_metrics
    }

}
