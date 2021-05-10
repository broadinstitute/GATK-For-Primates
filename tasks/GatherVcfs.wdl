version development



task HardFilteredVcfs {
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



task UnfilteredVcfs {
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