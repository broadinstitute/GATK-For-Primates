version 1.0

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
import "../tasks/ProcessVCFs.wdl" as processVCFs
import "../tasks/Utilities.wdl" as utilities

##########################################################################
## WORKFLOW DEFINITION
##########################################################################

workflow pipeGenotypes {

    input {
        File ref
        File ref_dict
        File ref_fai
        String gdb_groupName
        Array[File] gdb_vcfs
        Boolean pipe_unfiltered
    	Boolean pipe_unfiltered_sites_only
    	Boolean pipe_hard_filtered_SNPs
   		Boolean pipe_hard_filtered_SNPs_sites_only
    	Boolean pipe_hard_filtered_INDELs
    	Boolean pipe_hard_filtered_INDELs_sites_only
        String container
    }


	## If unfiltered variants are needed in any form
	if (pipe_unfiltered || pipe_unfiltered_sites_only) {

		## Get indexes
		scatter (each_gdb_vcf in gdb_vcfs) {
			String gdb_vcf_indexes = each_gdb_vcf + ".tbi"
		}

		## Gather unfiltered variants
        call processVCFs.gatherVCFs as gatherUnfiltered {
            input:
                input_vcfs = gdb_vcfs,
                input_vcfs_indexes = gdb_vcf_indexes,
                output_vcf_title = gdb_groupName + ".unfiltered_genotypes.",
                groupName = gdb_groupName,
                # Runtime
                container = container,
        }

        ## If unfiltered sites-only variants are needed
        if (pipe_unfiltered_sites_only) {

			## Scatter each file and make sites only
			scatter (each_gdb_vcf in gdb_vcfs) {

				String each_gdb_vcf_index = each_gdb_vcf + ".tbi"

	            call processVCFs.makeSitesOnly as makeSitesOnlyEachUnfiltered {
	                input:
	                    input_vcf = each_gdb_vcf,
	                    input_vcf_index = each_gdb_vcf_index,
	                    output_vcf_title = basename(each_gdb_vcf, ".vcf.gz") + ".unfiltered_genotypes_sites_only.",
	                    groupName = gdb_groupName,
	                    # Runtime
	                    container = container,
	            }
			}

			## Gather unfiltered variants
	        call processVCFs.gatherVCFs as makeSitesOnlyUnfiltered {
	            input:
	                input_vcfs = makeSitesOnlyEachUnfiltered.output_vcf,
	                input_vcfs_indexes = makeSitesOnlyEachUnfiltered.output_vcf_index,
	                output_vcf_title = gdb_groupName + ".unfiltered_genotypes.",
	                groupName = gdb_groupName,
	                # Runtime
	                container = container,
	        }

        }

    }



	## If hard-filtered SNPs or INDELs are needed in any form
	if (pipe_hard_filtered_SNPs || pipe_hard_filtered_SNPs_sites_only || pipe_hard_filtered_INDELs || pipe_hard_filtered_INDELs_sites_only) {

		scatter (gdb_vcf in gdb_vcfs) {

			String gdb_vcf_index = gdb_vcf + ".tbi"

			## Split multi-allelic rows
		    call processVCFs.splitMultiAllelics as splitMultiAllelics {
		        input:
		            ref = ref,
		            ref_dict = ref_dict,
		            ref_fai = ref_fai,
		            input_vcf = gdb_vcf,
		            input_vcf_index = gdb_vcf_index,
		            container = container,
		    }

		    ## Split into SNPs and INDELs
		    call processVCFs.splitVCFs as splitVCFs {
		        input:      
		            input_vcf = splitMultiAllelics.output_vcf,
		            input_vcf_index = splitMultiAllelics.output_vcf_index,
		            # Runtime
		            container = container,
		    }

		    ## If hard-filtered SNPs are needed in any form
			if (pipe_hard_filtered_SNPs || pipe_hard_filtered_SNPs_sites_only) {
				
			    call processVCFs.hardFilter as hardFilterSNPs {
			        input:
			        	ref = ref,
			        	ref_dict = ref_dict,
			        	ref_fai = ref_fai,
			            input_vcf = gdb_vcf,
			            input_vcf_index = gdb_vcf_index,
			            makeSitesOnly = pipe_hard_filtered_SNPs_sites_only,
			            type = "SNPs",
			            filters = "-filter 'QD < 2.0' --filter-name 'QD2' -filter 'QUAL < 30.0' --filter-name 'QUAL30' -filter 'SOR > 3.0' --filter-name 'SOR3' -filter 'FS > 60.0' --filter-name 'FS60' -filter 'MQ < 40.0' --filter-name 'MQ40' -filter 'MQRankSum < -12.5' --filter-name 'MQRankSum-12.5' -filter 'ReadPosRankSum < -8.0' --filter-name 'ReadPosRankSum-8'",
			            # Runtime
			            container = container,
			    }

			}

		    ## If hard-filtered INDELs are needed in any form
			if (pipe_hard_filtered_INDELs || pipe_hard_filtered_INDELs_sites_only) {

			    call processVCFs.hardFilter as hardFilterINDELs {
			        input:
			        	ref = ref,
			        	ref_dict = ref_dict,
			        	ref_fai = ref_fai,
			            input_vcf = gdb_vcf,
			            input_vcf_index = gdb_vcf_index,
						makeSitesOnly = pipe_hard_filtered_INDELs_sites_only,
						type = "INDELs",
			            filters = "-filter 'QD < 2.0' --filter-name 'QD2' -filter 'QUAL < 30.0' --filter-name 'QUAL30' -filter 'FS > 200.0' --filter-name 'FS200' -filter 'ReadPosRankSum < -20.0' --filter-name 'ReadPosRankSum-20'",
			            # Runtime
			            container = container,
			    }
			}

		}

		## Gather hard-filtered INDELs
		if (pipe_hard_filtered_INDELs || pipe_hard_filtered_INDELs_sites_only) {		
	        call processVCFs.gatherVCFs as gatherHFINDELs {
	            input:
	                input_vcfs = select_all(hardFilterINDELs.output_vcf),
	                input_vcfs_indexes = select_all(hardFilterINDELs.output_vcf_index),
	                output_vcf_title = gdb_groupName + ".hard_filtered_INDELs_only",
	                groupName = gdb_groupName,
	                # Runtime
	                container = container,
	        }

	        ## Gather hard-filtered sites only INDELs
			if (pipe_hard_filtered_INDELs_sites_only) {
	            call processVCFs.gatherVCFs as gatherHFINDELs_sites_only {
	                input:
	                    input_vcfs = select_all(hardFilterINDELs.output_sites_only_vcf),
	                    input_vcfs_indexes = select_all(hardFilterINDELs.output_sites_only_vcf_index),
	                    output_vcf_title = gdb_groupName + ".hard_filtered_INDELs_sites_only",
	                    groupName = gdb_groupName,
	                    # Runtime
	                    container = container,
	            }
			}
		}


		## Gather hard-filtered SNPs
		if (pipe_hard_filtered_SNPs || pipe_hard_filtered_SNPs_sites_only) {		
	        call processVCFs.gatherVCFs as gatherHFSNPs {
	            input:
	                input_vcfs = select_all(hardFilterSNPs.output_vcf),
	                input_vcfs_indexes = select_all(hardFilterSNPs.output_vcf_index),
	                output_vcf_title = gdb_groupName + ".hard_filtered_SNPs_only",
	                groupName = gdb_groupName,
	                # Runtime
	                container = container,
	        }

	        ## Gather hard-filtered sites only SNPs
			if (pipe_hard_filtered_SNPs_sites_only) {
	            call processVCFs.gatherVCFs as gatherHFSNPs_sites_only {
	                input:
	                    input_vcfs = select_all(hardFilterSNPs.output_sites_only_vcf),
	                    input_vcfs_indexes = select_all(hardFilterSNPs.output_sites_only_vcf_index),
	                    output_vcf_title = gdb_groupName + ".hard_filtered_SNPs_sites_only",
	                    groupName = gdb_groupName,
	                    # Runtime
	                    container = container,
	            }
			}
	    }
	}


    output {
    	String groupName = gdb_groupName
		File? unfiltered_vcf = gatherUnfiltered.output_vcf
		File? unfiltered_vcf_index = gatherUnfiltered.output_vcf_index
		File? unfiltered_sites_only_vcf = makeSitesOnlyUnfiltered.output_vcf
		File? unfiltered_sites_only_vcf_index = makeSitesOnlyUnfiltered.output_vcf_index
		File? hard_filtered_SNPs_vcf = gatherHFSNPs.output_vcf
		File? hard_filtered_SNPs_vcf_index = gatherHFSNPs.output_vcf_index
		File? hard_filtered_INDELs_vcf = gatherHFINDELs.output_vcf
		File? hard_filtered_INDELs_vcf_index = gatherHFINDELs.output_vcf_index
		File? hard_filtered_SNPs_sites_only_vcf = gatherHFSNPs_sites_only.output_vcf
		File? hard_filtered_SNPs_sites_only_vcf_index = gatherHFSNPs_sites_only.output_vcf_index
		File? hard_filtered_INDELs_sites_only_vcf = gatherHFINDELs_sites_only.output_vcf
		File? hard_filtered_INDELs_sites_only_vcf_index = gatherHFINDELs_sites_only.output_vcf_index
		String loc_unfiltered_vcf = select_first([gatherUnfiltered.output_vcf, "NULL"])
		String loc_unfiltered_vcf_index = select_first([gatherUnfiltered.output_vcf_index, "NULL"])
		String loc_unfiltered_sites_only_vcf = select_first([makeSitesOnlyUnfiltered.output_vcf, "NULL"])
		String loc_unfiltered_sites_only_vcf_index = select_first([makeSitesOnlyUnfiltered.output_vcf_index, "NULL"])
		String loc_hard_filtered_SNPs_vcf = select_first([gatherHFSNPs.output_vcf, "NULL"])
		String loc_hard_filtered_SNPs_vcf_index = select_first([gatherHFSNPs.output_vcf_index, "NULL"])
		String loc_hard_filtered_INDELs_vcf = select_first([gatherHFINDELs.output_vcf, "NULL"])
		String loc_hard_filtered_INDELs_vcf_index = select_first([gatherHFINDELs.output_vcf_index, "NULL"])
		String loc_hard_filtered_SNPs_sites_only_vcf = select_first([gatherHFSNPs_sites_only.output_vcf, "NULL"])
		String loc_hard_filtered_SNPs_sites_only_vcf_index = select_first([gatherHFSNPs_sites_only.output_vcf_index, "NULL"])
		String loc_hard_filtered_INDELs_sites_only_vcf = select_first([gatherHFINDELs_sites_only.output_vcf, "NULL"])
		String loc_hard_filtered_INDELs_sites_only_vcf_index = select_first([gatherHFINDELs_sites_only.output_vcf_index, "NULL"])
    }

}
