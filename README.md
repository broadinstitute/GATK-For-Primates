# GATK-for-Animals: GATK Best Practices for Variant Calling in Non-Human Animal Genomes

## IMPORTANT NOTE: THIS IS A PRE-ALPHA WORKFLOW CURRENTLY UNDER DEVELOPMENT!


### Purpose:

This workflow (gatk-for-animals-germline-snps-indels.wdl) facilitates [germline short variant discovery](https://gatk.broadinstitute.org/hc/en-us/articles/360035535932) from whole genomes of non-human animals, following proposed GATK Best Practices for Non-Human Animal Genomes. Capitalizing on new features in GATK4, the pipeline enables base recalibration and variant calling in the absence of 'gold standard' truth and training sets. Though principally designed for mammalian genomes, the pipeline may be adapted for other non-human animal species. Beyond animals, more appropriate GATK Best Practices might be available, e.g. [GATK for Microbes](https://github.com/broadinstitute/GATK-for-Microbes).

### Modes of operation:

The workflow operates in one of three modes, as specified in the mandatory input JSON file: initial, repeat or final.

#### 'Initial' mode
In initial mode, the workflow takes non-interleaved paired-end FASTQ reads, maps these to the reference genome, performs an initial round of variant calls, hard-filters those calls, and uses these to perform Base Quality Score Recalibration (BQSR). The output includes recalibrated BAM files for each individual, plus the necessary to tables and plots to evaluate if convergence has been reached.

#### 'Repeat' mode
In repeat mode, the workflow operates as in 'initial' mode, but starts with input BAM files instead of FASTQ reads. The input BAM files could either be a) BAM files that the user has previously mapped, and which require recalibration, or b) the output BAM files from 'initial' mode, which require further rounds of BQSR until covergence is reached.

#### 'Final' mode
In final mode, the workflow calls variants and then performs either Variant Quality Score Recalibration (VQSR) and/or hard-filtering of variants. The output comprises genomicsdb files (to which further samples can be added later, beyond this pipeline) plus final variant calls from downstream analysis.

### Truth and training sets

'Truth' sets for either SNPs and/or INDELs can be supplied in the input JSON. If both are supplied, the workflow will perform VQSR using hard-filtered variants as the training set and the input files at the truth sets.

If a 'truth' set is supplied for either SNPs or INDELs (but not both), then one will be subject to VQSR and the other to hard filtering. This is likely to be common in wildlife cases: SNP truth sets may be readily available (e.g. in the form of SNP microarray data) whereas INDEL truth sets are more difficult to come by.

If no truth sets are supplied, VQSR will be skipped and hard-filtering performed instead.


### Structure of input samples

In its simplest form, this pipeline can be run for multiple individuals comprising the same taxonomic unit (e.g. species). *NOTE: THIS DOESN'T WORK AT PRESENT*

However, a common-use scenario in wildlife comparative genomics is to call variants across individuals from multiple (but closely related) taxonomic units -- e.g. multiple species within a genus. In this case, the pipeline will call variants in the individuals of each taxonomic unit separately (e.g. each species), plus in all individuals combined (e.g. across the genus). This is done to identify rare alleles that might only persist within a single unit (e.g. in only one species) plus to identify low-frequency alleles (e.g. those across the genus). All these loci are then combined into a master callset, and each individual is re-genotyped at all of these loci. *NOTE: THIS WORKS AT PRESENT... HOWEVER, IT DOES NOT CALL ALL INDIVIDUALS TOGETHER (AS GENUS). I WILL IMPLEMENT THIS SOON.*

The pipeline will determine how best to operate based on the input JSON file. If individuals of only one unit (defined as 'group' in the JSON) are supplied, then the pipeline will perform in the former fashion. If multiple units are supplied, the pipeline will call separately for each unit and again for all samples together.

### Requirements

#### Input requirements/expectations:
- Paired-end reads in FASTQ format (for 'initial' mode) or mapped and de-duplicated BAM files with read group information (for 'repeat' and 'final' modes).
- A reference genome indexed using either bwa or bwa-mem2.
- An input JSON file listing the reference, the required metadata for each sample, and any optional parameters.
- Truth SNP and/or INDEL sets, if available.
- A list of 'scatters', either produced with ScatterIntervalsByNs (as described in the paper) or with the third-party Python tool, [ScaffoldStitcher](https://github.com/ameliahaj/ScaffoldStitcher)


#### Software version requirements :
- GATK 4.2.0.0
- Samtools 1.1.2
- Python 2.7
- Cromwell version support 
  - Testing in progress on v58

### IMPORTANT NOTES AND CAVEATS :
- The workflow requires the 'development' version of WDL.
- The workflow implements the GATK allele-specific filtering workflow that is currently in beta.
- At present, SNPs and INDELs are either recalibrated or filtered separately. I need to re-code this so that, if truth sets are supplied for SNPs and INDELs, they are recalibrated in series versus in sequence. Otherwise, [mixed sites will fail to be processed correctly in alelle-specific mode](https://gatk.broadinstitute.org/hc/en-us/articles/360035890551?id=9622).
- The --merge-contigs-into-num-partitions parameter is set to 0 (default), but can (and probably should) be configured to a higher value in the input JSON [as explained here](https://gatk.broadinstitute.org/hc/en-us/articles/360056138571-GDBI-usage-and-performance-guidelines). When preparing the final GenomicsDB during re-genotyping at the master callset of loci, the parameter is set to "1" (due to many small contigs, each comprising very little data).
- Re-genotyping at the master callset capitalizes on the aforementioned "--merge-contigs-into-num-partitions" parameter in ImportGenomicsDB, plus on the "--include-non-variant-sites" option in GenotypeGVCFs [as implemented in GATK 4.0.12.0](https://github.com/broadinstitute/gatk/releases/tag/4.0.12.0).
- Code for producing AnalyzeCovariates plots is currently commented out because the GATK4 docker does not have the required R libraries. Is there a docker image that has both R and the gsalib and ggplot2?
- No error checking/validation of input JSON options is implemented yet. Inputting contradictory options (e.g. mode = initial but with BAM inputs) may cause the workflow to fail.
- I had to remove "--merge-contigs-into-num-partitions" from GenomicsDBImport because it doesn't work as I expected -- it requires whole contigs versus intervals. I will re-visit this.
- Some tasks run on docker, some don't. I'll fix this.
- I haven't even tried this with Terra yet.
