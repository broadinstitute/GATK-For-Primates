version 1.0

## Copyright Broad Institute and Wisconsin National Primate Research Center,
## University of Wisconsin-Madison, 2021
## 
## Complete germline short variant discovery pipeline optimized for
## non-human primates, following proposed GATK Best Practices for
## Non-Human Animal Genomes.
##
## For requirements, expectations and outputs, please review the complete
## documentation at: https://github.com/broadinstitute/GATK-For-Primates
##
## Software version requirements :
## - Cromwell 65
## - bwa 0.7.15 (note: from GITC)
## - Samtools 1.11 (note: from GITC)
## - GATK 4.2.1.0 (note: GATK 4.1.8.0 is used in GITC)
## - Python 3.9.5
##
## Program versions can be changed by defining alternative containers.
## Runtime parameters are optimized for Terra (https://www.terra.bio/)
##
## LICENSING : 
## This script is released under the WDL source code license (BSD-3) (see LICENSE in 
## https://github.com/broadinstitute/wdl). However, the programs it calls may be
## subject to different licenses. Users are responsible for checking that they are
## authorized to run all programs before running this script. Please see the docker
## containers for detailed licensing information pertaining to the included programs.

import "structs/structs.wdl"
import "workflows/ValidateUserInputs.wdl" as preflight
import "workflows/UnmappedInputToAlignedBAM.wdl" as mapReads
import "tasks/Utilities.wdl" as utilities
import "tasks/HaplotypesToGenotypes.wdl" as bamToVcf
import "workflows/PipeGenotypes.wdl" as pipeGenotypes
import "tasks/ProcessVCFs.wdl" as processVCFs
import "workflows/VariantFiltration.wdl" as variantFiltration
import "workflows/VariantRecalibration.wdl" as variantRecalibration
import "workflows/VariantRecalibrationPartial.wdl" as partialVariantRecalibration
import "workflows/BaseRecalibration.wdl" as baseRecalibration
import "workflows/functions/CollectByKey.wdl" as collectByKey

##########################################################################
## WORKFLOW DEFINITION
##########################################################################

workflow GATKForPrimatesGermlineSNPsIndels_GATK4 {

    String pipeline_version = "alpha"

    ##########################################################################
    ## DEFINE WORFKLOW INPUTS
    ##########################################################################

    input {

        ## Collect workflow mode from input JSON
        String mode # options: initial / repeat / final
 
        ## Collect optional variables from input JSON
        Boolean validate_truth_sets = true # options: true / false; if false this will disable running ValidateVariants on truth sets in 'Final' mode
        Boolean flowcell_patterned = true # options: true / false; this influences pixel distance when marking duplicates
        Int? merge_contigs_into_num_partitions # options: optional parameter for GenomicsDBImport
        Boolean bwamem2 = false # options: true / false; indicating bwa (as bwa mem) or bwamem2 (as bwamem2 mem) ***-Coming-Soon-***
        Boolean cram_not_bam = true # options: true / false; if false this will disable use of CRAM instead of BAM format

        ## Collect reference files from input JSON
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

        ## Collect .tar.gz of packaged interval lists from 'initial' mode
        File? packaged_polymorphic_regions

        ## Collect truth sets from input JSON
        File? truth_set_SNPs # options: optional SNP truth set (sites-only VCF file) for VQSR (training set is produced via hard filtering)
        File? truth_set_SNPs_index # index for the above
        File? truth_set_INDELs # options: optional INDEL truth set (sites-only VCF file) for VQSR (training set is produced via hard filtering)
        File? truth_set_INDELs_index # index for the above
        
        ## Define arrays from input JSON; definitions are in the structs/structs.wdl file
        Array[sampleInfo]+ sampleList
        Array[scatterInfo]+ scatterList
        
        ## Define containers
        String container_gatk = "broadinstitute/gatk:4.2.1.0"
        String container_gitc = "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.5.7-2021-06-09_16-47-48Z"
        String container_python = "python:3.9.5"

        ## Define paths
        String path_to_gitc = "/usr/gitc/"
        String path_to_gitc_gatk = "/usr/gitc/gatk4/"

    }

    ##########################################################################
    ## Process input options and set default values
    ##########################################################################

    ## Assign correct index files to ref_idxs depending on version of bwa
    Array[File?]+ ref_idxs = if bwamem2 then [ref_amb, ref_ann, ref_pac, ref_0123, ref_bwt_2bit_64] else [ref_amb, ref_ann, ref_pac, ref_bwt, ref_sa]
    ## Set mode Booleans as true or false
    Boolean mode_is_initial = mode == "initial"
    Boolean mode_is_repeat = mode == "repeat"
    Boolean mode_is_final = mode == "final"
    ## Set Booleans determining which genotype files should be built and piped downstream
    Boolean pipe_unfiltered = if (defined(truth_set_INDELs)) || (defined(truth_set_SNPs)) && (mode_is_final) then true else false
    Boolean pipe_unfiltered_sites_only = if (pipe_unfiltered) then true else false
    Boolean pipe_hard_filtered_SNPs = if ((defined(truth_set_INDELs)) && (mode_is_final)) || ((!defined(truth_set_INDELs)) && (!defined(truth_set_SNPs)) && (mode_is_final)) then true else false
    Boolean pipe_hard_filtered_INDELs = if ((defined(truth_set_SNPs)) && (mode_is_final)) || ((!defined(truth_set_INDELs)) && (!defined(truth_set_SNPs)) && (mode_is_final)) then true else false
    Boolean pipe_hard_filtered_SNPs_sites_only = if ((defined(truth_set_INDELs)) && (!defined(truth_set_SNPs)) && (mode_is_final)) then false else true
    Boolean pipe_hard_filtered_INDELs_sites_only = if ((defined(truth_set_SNPs)) && (!defined(truth_set_INDELs)) && (mode_is_final)) then false else true

    ##########################################################################
    ## Validate user input files and variables
    ##########################################################################

    call preflight.validateUserInputs as validateUserInputs {
        input:
            ref = ref,
            ref_dict = ref_dict,
            ref_fai = ref_fai,
            scatterList = scatterList,
            sampleList = sampleList,
            truth_set_SNPs = truth_set_SNPs,
            truth_set_SNPs_index = truth_set_SNPs_index,
            truth_set_INDELs = truth_set_INDELs,
            truth_set_INDELs_index = truth_set_INDELs_index,
            mode = mode,
            packaged_polymorphic_regions = packaged_polymorphic_regions,
            ## Runtime
            container_gatk = container_gatk,
            container_python = container_python,
    }

    Array[polymorphicRegionsInfo]+ polymorphicRegionsList = read_json(validateUserInputs.polymorphic_regions_json)

    ####################################################################
    ## IF MODE == INITIAL: Map either paired FASTQ or uBAM to BAM
    ####################################################################
    
    if (mode_is_initial) {

        call mapReads.unmappedInputToAlignedBAM as mapReads {
            input:
                ref = ref,
                ref_dict = ref_dict,
                ref_fai = ref_fai,
                ref_idxs = ref_idxs,
                sampleList = sampleList,
                flowcell_patterned = flowcell_patterned,
                bwamem2 = bwamem2,
                cram_not_bam = cram_not_bam,
                ## Runtime
                container_gatk = container_gatk,
                container_gitc = container_gitc,
                path_to_gitc = path_to_gitc,
                path_to_gitc_gatk = path_to_gitc_gatk,
                container_python = container_python,
        }

    }

    ####################################################################
    ## Generate JSON of BAM files to use downstream
    ####################################################################

    call utilities.generateBamJSON as bamJSON {
        input:
            sampleGroups = select_first([mapReads.sampleGroups,validateUserInputs.groupNames]),
            sampleNames = select_first([mapReads.sampleNames,validateUserInputs.sampleNames]),
            bams = select_first([mapReads.bams,validateUserInputs.bams]),
            bam_indexes = select_first([mapReads.bam_indexes,validateUserInputs.bam_indexes]),
            ## Runtime
            container = container_python,
    }

    Array[bamInfo]+ bamList = read_json(bamJSON.file)

    ####################################################################
    ## IF MODE == INITIAL: Call variants across user-defined scatters
    ####################################################################

    if (mode_is_initial) { 

        ## Scatter over each user-defined scatter interval
        scatter (scttr in scatterList) {

            ## Scatter over each BAM file
            scatter (bamFile in bamList) {

                ## Perform initial round of variant calling
                call bamToVcf.haplotypeCallerVCF as haplotypeCallerVCF {
                    input:
                        ref = ref,
                        ref_dict = ref_dict,
                        ref_fai = ref_fai,
                        input_bam = bamFile.bam,
                        input_bam_index = bamFile.bam_index,
                        sampleGroup = bamFile.sampleGroup,
                        sampleName = bamFile.sampleName,
                        scatterName = scttr.name,
                        scatterIntervals = sub(scttr.intervals, " ", " -L "),
                        # Runtime
                        container = container_gatk,
                }

                ## Split, hard filter and make sites only high-confidence SNPs and INDELs only
                call variantFiltration.variantFiltration as filterInitialVariants {
                    input:
                        ref = ref,
                        ref_dict = ref_dict,
                        ref_fai = ref_fai,
                        input_vcf = haplotypeCallerVCF.output_vcf,
                        input_vcf_index = haplotypeCallerVCF.output_vcf_index,
                        # Runtime
                        container = container_gatk,
                }

                ## Make unfiltered haplotypes to sites only
                call processVCFs.makeSitesOnly as makeSitesOnlyInitialVariants {
                    input:
                        groupName = "unfiltered",
                        input_vcf = haplotypeCallerVCF.output_vcf,
                        input_vcf_index = haplotypeCallerVCF.output_vcf_index,
                        # Runtime
                        container = container_gatk,
                } 

            }

            ## Produce and pre-process interval lists for high-confidence regions per scatter
            call processVCFs.polymorphicSitesToRegions as polymorphicSitesToRegions {
                input:
                    ref = ref,
                    ref_dict = ref_dict,
                    ref_fai = ref_fai,
                    scatterName = scttr.name,
                    input_vcf_array = select_all(makeSitesOnlyInitialVariants.output_vcf),
                    input_vcf_index_array = select_all(makeSitesOnlyInitialVariants.output_vcf_index),
                    # Runtime
                    container = container_gatk,
            }

        }

        ## Generate JSON of polymorphic regions across scatters
        call utilities.packagePolymorphicRegions as packagePolymorphicRegions {
            input:
                scatterNames = select_all(polymorphicSitesToRegions.output_scatterName),
                intervalLists = select_all(polymorphicSitesToRegions.output_intervalList),
                intervalListFilenames = select_all(polymorphicSitesToRegions.output_intervalList_filename),
                # Runtime
                container = container_python,
        }

    }


    ####################################################################
    ## IF MODE == REPEAT: Call variants across polymorphic regions
    ####################################################################

    if (mode_is_repeat) { 

        ## Scatter over each user-defined scatter interval
        scatter (scttrPoly in polymorphicRegionsList) {

            ## Fetch the correct interval list file for this scatter from the glob
            scatter (listLocation in validateUserInputs.polymorphicRegionsIntervalLists) {
                if (basename(listLocation) == scttrPoly.intervalList) {
                    String? listToUse = listLocation
                }
            }

            ## Scatter over each BAM file
            scatter (bamFile in bamList) {

                ## Perform initial round of variant calling
                call bamToVcf.haplotypeCallerVCF as repeatHaplotypeCallerVCF {
                    input:
                        ref = ref,
                        ref_dict = ref_dict,
                        ref_fai = ref_fai,
                        input_bam = bamFile.bam,
                        input_bam_index = bamFile.bam_index,
                        sampleGroup = bamFile.sampleGroup,
                        sampleName = bamFile.sampleName,
                        scatterName = scttrPoly.scatterName,
                        scatterIntervalList = select_first(listToUse),
                        # Runtime
                        container = container_gatk,
                }

                ## Split, hard filter and make sites only high-confidence SNPs and INDELs
                call variantFiltration.variantFiltration as filterRepeatVariants {
                    input:
                        ref = ref,
                        ref_dict = ref_dict,
                        ref_fai = ref_fai,
                        input_vcf = repeatHaplotypeCallerVCF.output_vcf,
                        input_vcf_index = repeatHaplotypeCallerVCF.output_vcf_index,
                        # Runtime
                        container = container_gatk,
                }

            }

        }

    }

    ####################################################################
    ## IF MODE == INITIAL/REPEAT: Recalibrate BAM files
    ####################################################################    

    if (mode_is_initial || mode_is_repeat) {

            Array[Array[File]] BQSR_SNP_sites_array = select_first([filterInitialVariants.SNP_sites_vcf,filterRepeatVariants.SNP_sites_vcf])
            Array[Array[File]] BQSR_SNP_sites_index_array = select_first([filterInitialVariants.SNP_sites_vcf_index,filterRepeatVariants.SNP_sites_vcf_index])
            Array[Array[File]] BQSR_INDEL_sites_array = select_first([filterInitialVariants.INDEL_sites_vcf,filterRepeatVariants.INDEL_sites_vcf])
            Array[Array[File]] BQSR_INDEL_sites_index_array = select_first([filterInitialVariants.INDEL_sites_vcf_index,filterRepeatVariants.INDEL_sites_vcf_index])

        ## Merge all filtered variants from all samples within scatter
        call processVCFs.mergeManyVCFs as mergeSitesForBQSR {
            input:
                ref_dict = ref_dict,
                input_array_1_vcf = flatten(BQSR_SNP_sites_array),
                input_array_1_vcf_index = flatten(BQSR_SNP_sites_index_array),
                input_array_2_vcf = flatten(BQSR_INDEL_sites_array),
                input_array_2_vcf_index = flatten(BQSR_INDEL_sites_index_array),
                output_vcf_title = "high_confidence_sites_for_BQSR",
                # Runtime
                container = container_gatk
        }

        ## Scatter over each input BAM file to recalibrate
        scatter (bamFile in bamList) {

            ## Perform and apply base recalibration
            call baseRecalibration.BQSR as BQSR {
                input:
                    ref = ref,
                    ref_dict = ref_dict,
                    ref_fai = ref_fai,
                    input_bam = bamFile.bam,
                    input_bam_index = bamFile.bam_index,
                    cram_not_bam = cram_not_bam,
                    high_confidence_sites_vcf = mergeSitesForBQSR.output_merged_vcf,
                    high_confidence_sites_vcf_index = mergeSitesForBQSR.output_merged_vcf_index,
                    sampleName = bamFile.sampleName,
                    container_gatk = container_gatk,
                    container_gitc = container_gitc,
                    path_to_gitc_gatk = path_to_gitc_gatk,
            }

        }

    }

    ####################################################################
    ## IF MODE == FINAL: 
    ####################################################################
    ## Call haplotypes and perform joint genotyping per taxon
    ####################################################################

    if (mode_is_final) { 

        ## Scatter over each scatter and polymorphic region
        scatter (scttr in polymorphicRegionsList) {

            ## Fetch the correct interval list file for this scatter from the glob
            scatter (listLocation in validateUserInputs.polymorphicRegionsIntervalLists) {
                if (basename(listLocation) == scttr.intervalList) {
                    String? listToUse_final = listLocation
                }
            }

            ## Scatter over each BAM file
            scatter (bamFile in bamList) {

                ## Call haplotypes (i.e. per BAM, per scatter)
                call bamToVcf.haplotypeCallerGVCF as haplotypeCallerGVCF {
                    input:
                        ref = ref,
                        ref_dict = ref_dict,
                        ref_fai = ref_fai,
                        input_bam = bamFile.bam,
                        input_bam_index = bamFile.bam_index,
                        sampleName = bamFile.sampleName,
                        sampleGroup = bamFile.sampleGroup,
                        scatterName = scttr.scatterName,
                        intervalList = select_first(listToUse_final),
                        # Runtime
                        container = container_gatk,
                }

            }

            ## Generate haplotype JSON per scatter, adding 'Cohort' group if multiple taxa
            call utilities.generateHaplotypeJSON as haplotypeJSON {
                input:
                    sampleGroups = haplotypeCallerGVCF.output_sampleGroup,
                    gvcfs = haplotypeCallerGVCF.output_gvcf,
                    multiple_taxonomic_groups = validateUserInputs.multiple_taxonomic_groups,
                    # Runtime
                    container = container_python,
            }

            ## Read haplotypeJSON as a haplotypeInfo struct and coerce into arrays
            Array[haplotypeInfo]+ haplotypeList = read_json(haplotypeJSON.file)
            scatter (haplos in haplotypeList) {
                String haplotype_groupNames = haplos.groupName
                String haplotype_gvcfs = haplos.gvcf
            }

            ## Scatter over each group, collecting haplotype gVCF files for each
            ## WDL 1.1: scatter (pair in as_pairs(collect_by_key(zip(select_all(haplotype_groupNames),select_all(haplotype_gvcfs))))) {
            call collectByKey.collectByKey as collectHaplotypesByGroup {
                input:
                    allGroups = validateUserInputs.groupNames,
                    groups = select_all(haplotype_groupNames),
                    members = select_all(haplotype_gvcfs),
            }

            scatter (pair in collectHaplotypesByGroup.collected) {

                scatter (each_haplotype_gvcf in pair.right) {
                    String haplotype_gvcf_indexes = each_haplotype_gvcf + ".tbi"
                }
                
                ## Import haplotypes to GenomicsDBs
                call bamToVcf.genomicsDBImport as genomicsDBImport {
                    input:
                        input_gvcfs = pair.right,
                        input_gvcfs_indexes = haplotype_gvcf_indexes,
                        sampleGroup = pair.left,
                        scatterName = scttr.scatterName,
                        intervalList = select_first(listToUse_final),
                        merge_contigs_into_num_partitions = merge_contigs_into_num_partitions,
                        # Runtime
                        container = container_gatk,
                }

                ## Genotype each database across polymorphic regions
                call bamToVcf.genotypeGenomicsDB as genotypeGenomicsDB {
                    input:
                        ref = ref,
                        ref_dict = ref_dict,
                        ref_fai = ref_fai,
                        input_genomicsdb = genomicsDBImport.output_genomicsdb,
                        input_genomicsdb_title = genomicsDBImport.output_genomicsdb_title,
                        output_vcf_title = "variant_calls_" + pair.left + "_" + scttr.scatterName + "_genotyped",
                        groupName = pair.left,
                        scatterName = scttr.scatterName,
                        intervalList = select_first(listToUse_final),
                        regenotype = false,
                        # Runtime
                        container = container_gatk,
                }

            }

        }

        ####################################################################
        ## Produce and pipe the variants that are needed downstream
        ####################################################################

        ## Scatter over each group, collecting genotype VCF files for each
        ## WDL 1.1: scatter (pair in as_pairs(collect_by_key(zip(select_all(flatten(genotypeGenomicsDB.output_groupName)),select_all(flatten(genotypeGenomicsDB.output_genotypes_vcf)))))) {
        call collectByKey.collectByKey as collectGenotypesByGroup {
            input:
                allGroups = validateUserInputs.groupNames,
                groups = select_all(flatten(genotypeGenomicsDB.output_groupName)),
                members = select_all(flatten(genotypeGenomicsDB.output_genotypes_vcf)),
        }
        Array[Pair[String,Array[String]]] collectedGenotypesByGroup = collectGenotypesByGroup.collected
    
        scatter (pair in collectedGenotypesByGroup) {

            call pipeGenotypes.pipeGenotypes as pipe {
                input:
                    ref = ref,
                    ref_dict = ref_dict,
                    ref_fai = ref_fai,
                    gdb_groupName = pair.left,
                    gdb_vcfs = pair.right,
                    pipe_unfiltered = pipe_unfiltered,
                    pipe_unfiltered_sites_only = pipe_unfiltered_sites_only,
                    pipe_hard_filtered_SNPs = pipe_hard_filtered_SNPs,
                    pipe_hard_filtered_INDELs = pipe_hard_filtered_INDELs,
                    pipe_hard_filtered_SNPs_sites_only = pipe_hard_filtered_SNPs_sites_only,
                    pipe_hard_filtered_INDELs_sites_only = pipe_hard_filtered_INDELs_sites_only,
                    # Runtime
                    container = container_gatk,
            }

        }

        call utilities.generateGenotypeJSON as genotypeJSON {
            input:
                groupName = pipe.groupName,
                loc_unfiltered_vcf = pipe.loc_unfiltered_vcf,
                loc_unfiltered_vcf_index  = pipe.loc_unfiltered_vcf_index,
                loc_unfiltered_sites_only_vcf = pipe.loc_unfiltered_sites_only_vcf,
                loc_unfiltered_sites_only_vcf_index = pipe.loc_unfiltered_sites_only_vcf_index,
                loc_hard_filtered_SNPs_vcf = pipe.loc_hard_filtered_SNPs_vcf,
                loc_hard_filtered_SNPs_vcf_index = pipe.loc_hard_filtered_SNPs_vcf_index,
                loc_hard_filtered_INDELs_vcf = pipe.loc_hard_filtered_INDELs_vcf,
                loc_hard_filtered_INDELs_vcf_index = pipe.loc_hard_filtered_INDELs_vcf_index,
                loc_hard_filtered_SNPs_sites_only_vcf = pipe.loc_hard_filtered_SNPs_sites_only_vcf,
                loc_hard_filtered_SNPs_sites_only_vcf_index = pipe.loc_hard_filtered_SNPs_sites_only_vcf_index,
                loc_hard_filtered_INDELs_sites_only_vcf = pipe.loc_hard_filtered_INDELs_sites_only_vcf,
                loc_hard_filtered_INDELs_sites_only_vcf_index = pipe.loc_hard_filtered_INDELs_sites_only_vcf_index,
                # Runtime
                container = container_python,
        }

        Array[genotypeInfo]+ genotypeList = read_json(genotypeJSON.file)

        ####################################################################
        ## No truth sets available? Generate hard-filtered callset...
        ####################################################################

        if ((!defined(truth_set_INDELs)) && (!defined(truth_set_SNPs))) {

            call processVCFs.mergeManyVCFs as generateHardFilteredCallset {
                input:
                    input_array_1_vcf = select_all(pipe.hard_filtered_SNPs_sites_only_vcf),
                    input_array_1_vcf_index = select_all(pipe.hard_filtered_SNPs_sites_only_vcf_index),
                    input_array_2_vcf = select_all(pipe.hard_filtered_INDELs_sites_only_vcf),
                    input_array_2_vcf_index = select_all(pipe.hard_filtered_INDELs_sites_only_vcf_index),
                    ref_dict = ref_dict,
                    output_vcf_title = "hardFilteredCallset",
                    # Runtime
                    container = container_gatk,
            }
            
        }

        ####################################################################
        ## One or both truth sets available? Proceed with VQSR...
        ####################################################################

        ## If mode is final and at least one truth set is available
        if ((defined(truth_set_INDELs)) || (defined(truth_set_SNPs))) {

            ## Scatter over the available genotype sets per group
            scatter (genotypeFile in genotypeList) {

                ## If both truth sets are available, recalibrate each in series
                if (defined(truth_set_SNPs) && defined(truth_set_INDELs)) {
                    call variantRecalibration.VQSR as VQSR {
                        input:
                            ref = ref,
                            ref_dict = ref_dict,
                            ref_fai = ref_fai,
                            groupName = genotypeFile.groupName,
                            unfiltered_vcf = genotypeFile.unfiltered_vcf,
                            unfiltered_vcf_index = genotypeFile.unfiltered_vcf_index,
                            unfiltered_vcf_sites_only = genotypeFile.unfiltered_sites_only_vcf,
                            unfiltered_vcf_sites_only_index = genotypeFile.unfiltered_sites_only_vcf_index,
                            truth_set_INDELs = truth_set_INDELs,
                            truth_set_INDELs_index = truth_set_INDELs_index,
                            truth_set_SNPs = truth_set_SNPs,
                            truth_set_SNPs_index = truth_set_SNPs_index,
                            training_set_INDELs = genotypeFile.hard_filtered_INDELs_sites_only_vcf,
                            training_set_INDELs_index = genotypeFile.hard_filtered_INDELs_sites_only_vcf_index,
                            training_set_SNPs = genotypeFile.hard_filtered_SNPs_sites_only_vcf,
                            training_set_SNPs_index = genotypeFile.hard_filtered_SNPs_sites_only_vcf_index,
                            container = container_gatk,
                        }
                }

                ## If only one is available, perform partial recalibration combined with hard filtering
                if (!defined(truth_set_SNPs) || !defined(truth_set_INDELs)) {
                    call partialVariantRecalibration.partialVQSR as partialVQSR {
                        input:
                            ref = ref,
                            ref_dict = ref_dict,
                            ref_fai = ref_fai,
                            groupName = genotypeFile.groupName,
                            unfiltered_vcf = genotypeFile.unfiltered_vcf,
                            unfiltered_vcf_index = genotypeFile.unfiltered_vcf_index,
                            unfiltered_vcf_sites_only = genotypeFile.unfiltered_sites_only_vcf,
                            unfiltered_vcf_sites_only_index = genotypeFile.unfiltered_sites_only_vcf_index,
                            type = if defined(truth_set_SNPs) && !defined(truth_set_INDELs) then "SNP" else "INDEL",
                            truth_set = if defined(truth_set_SNPs) && !defined(truth_set_INDELs) then truth_set_SNPs else truth_set_INDELs,
                            truth_set_index = if defined(truth_set_SNPs) && !defined(truth_set_INDELs) then truth_set_SNPs_index else truth_set_INDELs_index,
                            training_set = if defined(truth_set_SNPs) && !defined(truth_set_INDELs) then genotypeFile.hard_filtered_SNPs_sites_only_vcf else genotypeFile.hard_filtered_INDELs_sites_only_vcf,
                            training_set_index = if defined(truth_set_SNPs) && !defined(truth_set_INDELs) then genotypeFile.hard_filtered_SNPs_sites_only_vcf_index else genotypeFile.hard_filtered_INDELs_sites_only_vcf_index,
                            container = container_gatk,
                        }
                }

            } ## Stop scattering over genotype sets per group
            
            ## If only one is available, merge the partially VQSR-filtered variants from all groups with the opposite hard-filtered variants from all groups
            if (!defined(truth_set_SNPs) || !defined(truth_set_INDELs)) {

                ## VQSR'd INDELs, hard-filtered SNPs
                if (!defined(truth_set_SNPs) && defined(truth_set_INDELs)) {
                    call processVCFs.mergeManyVCFs as mergeINDELsFromVQSRwithHardFilteredSNPs {
                        input:
                            input_array_1_vcf = select_all(partialVQSR.output_vcf),
                            input_array_1_vcf_index = select_all(partialVQSR.output_vcf_index),
                            input_array_2_vcf = select_all(pipe.hard_filtered_SNPs_sites_only_vcf),
                            input_array_2_vcf_index = select_all(pipe.hard_filtered_SNPs_sites_only_vcf_index),
                            ref_dict = ref_dict,
                            output_vcf_title = "mergeINDELsFromVQSRwithHardFilteredSNPs",
                            # Runtime
                            container = container_gatk,
                    }              
                }
                ## VQSR'd SNPs, hard-filtered INDELs
                if (defined(truth_set_SNPs) && !defined(truth_set_INDELs)) {
                    call processVCFs.mergeManyVCFs as mergeSNPsFromVQSRwithHardFilteredINDELs {
                        input:
                            input_array_1_vcf = select_all(partialVQSR.output_vcf),
                            input_array_1_vcf_index = select_all(partialVQSR.output_vcf_index),
                            input_array_2_vcf = select_all(pipe.hard_filtered_INDELs_sites_only_vcf),
                            input_array_2_vcf_index = select_all(pipe.hard_filtered_INDELs_sites_only_vcf_index),
                            ref_dict = ref_dict,
                            output_vcf_title = "mergeSNPsFromVQSRwithHardFilteredINDELs",
                            # Runtime
                            container = container_gatk,
                    }
                }

            }

            ## If both are available, merge the fully VQSR-filtered variants from all groups
            if (defined(truth_set_SNPs) && defined(truth_set_INDELs)) {

                # Merge variants into final loci callset
                call processVCFs.mergeVCFs as mergeAllFromVQSR {
                    input:
                        input_vcf_array = select_all(VQSR.output_vcf),
                        input_vcf_index_array = select_all(VQSR.output_vcf_index),
                        ref_dict = ref_dict,
                        output_vcf_title = "mergeAllFromVQSR",
                        # Runtime
                        container = container_gatk,
                }

            }

        }

        ####################################################################
        ## Re-genotype all samples at all loci
        ####################################################################

        ## Collect groupName as pair.left and genomicsDB as pair.right
        call collectByKey.collectByKey as collectGenomicsDBsByGroup {
            input:
                allGroups = validateUserInputs.groupNames,
                groups = select_all(flatten(genomicsDBImport.output_groupName)),
                members = select_all(flatten(genomicsDBImport.output_genomicsdb)),
        }

        if (validateUserInputs.multiple_taxonomic_groups) {
            String group_to_collect_mult = "Cohort"
        }

        if (!validateUserInputs.multiple_taxonomic_groups) {
            String group_to_collect_non_mult = select_first(flatten(genomicsDBImport.output_groupName))
        }

        String group_to_collect = select_first([group_to_collect_mult,group_to_collect_non_mult])

        scatter (pair in collectGenomicsDBsByGroup.collected) {
            if (pair.left == group_to_collect) {
                scatter (each_gdb in pair.right) {
                    # Genotype each database across polymorphic regions
                    call bamToVcf.genotypeGenomicsDB as genotypeGenomicsDBFinal {
                        input:
                            ref = ref,
                            ref_dict = ref_dict,
                            ref_fai = ref_fai,
                            input_genomicsdb = each_gdb,
                            input_genomicsdb_title = "FinalCallset",
                            output_vcf_title = "final_calls_" + basename(each_gdb),
                            groupName = group_to_collect,
                            scatterName = "FinalCallset",
                            intervalList = select_first([generateHardFilteredCallset.output_merged_vcf,mergeINDELsFromVQSRwithHardFilteredSNPs.output_merged_vcf,mergeSNPsFromVQSRwithHardFilteredINDELs.output_merged_vcf,mergeAllFromVQSR.output_merged_vcf]),
                            finalCallsetIndex = select_first([generateHardFilteredCallset.output_merged_vcf_index,mergeINDELsFromVQSRwithHardFilteredSNPs.output_merged_vcf_index,mergeSNPsFromVQSRwithHardFilteredINDELs.output_merged_vcf_index,mergeAllFromVQSR.output_merged_vcf_index]),
                            regenotype = true,
                            # Runtime
                            container = container_gatk,
                    }
                }
            }
        }

        call processVCFs.mergeVCFs as mergeFinalGenotypes {
            input:
                input_vcf_array = flatten(select_all(genotypeGenomicsDBFinal.output_genotypes_vcf)),
                input_vcf_index_array = flatten(select_all(genotypeGenomicsDBFinal.output_genotypes_vcf_index)),
                ref_dict = ref_dict,                
                output_vcf_title = "FinalGenotypes",
                # Runtime
                container = container_gatk,
        }




    }

    output {
        ## Outputs from initial mode:
        File? polymorphic_sites_JSON = packagePolymorphicRegions.file
        File? polymorphic_sites_tar = packagePolymorphicRegions.file_tar
        ## Outputs from initial/repeat modes:
        File? high_confidence_sites_BQSR_vcf = mergeSitesForBQSR.output_merged_vcf
        File? high_confidence_sites_BQSR_vcf_index = mergeSitesForBQSR.output_merged_vcf_index
        Array[File]? recalibrated_bam = BQSR.recalibrated_bam
        Array[File]? recalibrated_bam_index = BQSR.recalibrated_bam_index
        Array[File]? table_before = BQSR.table_before
        Array[File]? table_after = BQSR.table_after
        Array[File]? plots = BQSR.plots
        ## Outputs from final mode
        #File? finalGenotypes = mergeFinalGenotypes.output_merged_vcf
        #File? finalGenotypesIndex = mergeFinalGenotypes.output_merged_vcf_index
        #File? finalCallset = gatherFinalCallset.output_vcf
        #File? finalCallsetIndex = gatherFinalCallset.output_vcf_index
    }
    
}
