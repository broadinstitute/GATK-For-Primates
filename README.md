# GATK-for-Primates: GATK Best Practices for Variant Calling in Non-Human Primate Genomes (alpha)

## This repo is NOT ready for use yet and is still under development!!

A reproducible pipeline for germline SNP and Indel variant calling in non-human primate whole-genome re-sequencing data.

This pipeline is in the **alpha** stage of development. Please do not use this pipeline until we have completed initial benchmarking. Once we see promising results with sensitivity and precision, we'll seek the community's feedback to make additional improvements.

**Run on [Cromwell](https://cromwell.readthedocs.io/en/stable/):** `cromwell run gatk-for-primates-germline-snps-indels.wdl -i inputs.json`

**Run on [Terra](https://terra.bio/):** See the `GATK_For_Primates` workspace (coming soon: stay tuned!).

## Summary:

The [Genome Analysis Toolkit (GATK)](https://gatk.broadinstitute.org/) remains the premier software package for variant discovery and genotyping from next-generation resequencing reads. Though originally developed for human genetics, GATK has evolved to handle genome data from any organism with any level of ploidy. Nonetheless, the existing [‘best practices’](https://currentprotocols.onlinelibrary.wiley.com/doi/abs/10.1002/0471250953.bi1110s43) and [published workflows](https://gatk.broadinstitute.org/hc/en-us/sections/360007226651-Best-Practices-Workflows) are focused solely on processing human data. These cannot accommodate the unique challenges encountered in data from other taxa. Here, we present a 'best practices' [WDL pipeline](https://gatk.broadinstitute.org/hc/en-us/articles/360035889771-Pipelining-GATK-with-WDL-and-Cromwell) for [germline short variant discovery](https://gatk.broadinstitute.org/hc/en-us/articles/360035535932) in non-human primate whole genome re-sequencing data.

* This pipeline fully automates [Base Quality Score Recalibration (BQSR)](https://gatk.broadinstitute.org/hc/en-us/articles/360035890531-Base-Quality-Score-Recalibration-BQSR-) by applying hard-filtered variants from initial rounds of variant calling. This replaces the need for 'gold standard' [sets of known-variant sites](https://gatk.broadinstitute.org/hc/en-us/articles/360035890811-Resource-bundle) that only derive from humans.

* [Hard-filtering is also used](https://gatk.broadinstitute.org/hc/en-us/articles/360035531112--How-to-Filter-variants-either-with-VQSR-or-by-hard-filtering) in place of [Variant Quality Score Recalibration (VQSR)](https://gatk.broadinstitute.org/hc/en-us/articles/360035531612-Variant-Quality-Score-Recalibration-VQSR-), though partial or full VQSR will proceed if SNP and/or Indel truthing sets are provided.

* Samples from multiple (but closely related) taxonomic units (_e.g_. species within a genus, or subspecies within a species) can be co-analysed together in this pipeline. This allows for downstream analyses that, while irrelevant to human studies, transpire to be common-use scenarios in wildlife comparative genomics (_e.g._ assessing hybridization, admixture, speciation, _etc._). This approach also facilitates the discovery and consistent genotyping of rare and low-frequency alleles across taxonomic units, which may offer improved insights into complex evolutionary histories and traits.

* To ensure run-time efficiency, we take a [‘scatter-gather’ approach to parallelization](https://gatk.broadinstitute.org/hc/en-us/articles/360035532012-Parallelism-Multithreading-Scatter-Gather), which can be easily adapted for use on single computers or any high-throughput platform. The pipeline will also be available as a workspace on [Terra](https://www.terra.bio/) shortly.

Though presently under testing in primates, the pipeline may be adapted for other non-human animal species, for which new GATK Best Practices are currently in development. Beyond animals, more appropriate GATK Best Practices might be available, e.g. via the [GATK for Microbes](https://github.com/broadinstitute/GATK-for-Microbes) pipeline.

## Modes of operation:

The pipeline operates in one of three modes, to be specified in the input JSON file: `initial`, `repeat` or `final`. These comprise the major stages of a typical workflow, _i.e._ 1) mapping and BQSR, 2) repeat rounds of BQSR and 3) genotyping and VQSR and/or variant filtration. First-time users may want to run these sequentially, starting with raw FASTQ files and ending with a VCF. Returning and advanced users can simply change the `mode` to begin at a later stage.

#### Initial mode
In `initial` mode, the pipeline takes paired-end FASTQ reads or unmapped BAM files, **1)** maps the reads to the reference genome, **2)** performs an initial round of variant calling using HaplotypeCaller, **3)** hard-filters those calls using [standard parameters](https://gatk.broadinstitute.org/hc/en-us/articles/360035531112), and **4)** uses the passing loci to perform BQSR. The unfiltered variant calls are then processed into [Picard-style interval lists](https://gatk.broadinstitute.org/hc/en-us/articles/360035531852-Intervals-and-interval-lists) (which are packaged as `polymorphicRegions.tar.gz`) and outputted for use in `repeat` and `final` modes. The pipeline also outputs recalibrated BAM files for each individual, plus [the tables and plots required to evaluate convergence](https://gatk.broadinstitute.org/hc/en-us/articles/360035890531-Base-Quality-Score-Recalibration-BQSR-).

#### Repeat mode
`Repeat` mode is only used when there are no 'gold standard' [sets of known-variant sites](https://gatk.broadinstitute.org/hc/en-us/articles/360035890811-Resource-bundle). The mode should be run repeatedly, essentially 'bootstrapping' the BQSR step, until convergence (see examples [here](https://gatk.broadinstitute.org/hc/en-us/articles/360035890531-Base-Quality-Score-Recalibration-BQSR-), though in practice, this doesn't usually require more than two rounds). `Repeat` mode operates largely as in `initial` mode, except the user provides **(a)** the recalibrated BAM files (as outputted either from `initial` mode, or from the most recent run in `repeat` mode) and **(b)** the `polymorphicRegions.tar.gz` file (outputted from `initial` mode) as inputs. Instead of calling variants across the entirety of each BAM file, `repeat` mode restricts variant calling to the known (or at least suspected!) polymorphic regions, as indicated by the interval lists in the tarball. This is faster and far less computationally expensive than re-calling haplotypes across the entire genome. Previously mapped BAM files (_i.e._ from other sources) can also be provided at this stage, in lieu of unmapped data.

![image](https://user-images.githubusercontent.com/58449659/127539984-68e6b9ea-5ef0-4a49-baba-690f49c59e45.png)

#### Final mode
In `final` mode, the outputs from `initial` or `repeat` mode are used as inputs, _i.e._ fully recalibrated BAM files and the ``polymorphicRegions.tar.gz`` file comprising the interval lists of polymorphic regions. Variants are then called to intermediate GVCF files as follows:

- If samples from multiple taxonomic groups are present, each group's samples are first joint-genotyped separately, then again all together across the entire cohort.
- If 'gold standard' SNP and Indel sets are provided, these are used as 'truth' sets in VQSR, alongside hard-filtered 'training' sets generated by the pipeline. Recalibration is subsequently performed in series (Indels, then SNPs) across all unfiltered variant calls.
- If only one 'gold standard' set is provided (_e.g._ SNPs), recalibration is first performed on all unfiltered variant calls using that truth set; the results are then merged with the hard-filtered calls from the other variant type (_i.e._ Indels, in this example).

At the conclusion of `final` mode, all samples are re-genotyped at the passing loci. The final output includes [GenomicsDB](https://gatk.broadinstitute.org/hc/en-us/articles/360035891051-GenomicsDB) of all samples (to which further samples can be added later, beyond this pipeline) plus the final genotypes (in gzipped VCF format) for downstream applications.

![image](https://user-images.githubusercontent.com/58449659/127561880-41a375d9-abc7-4dcf-b216-7120f516235e.png)

## Requirements

### Software version requirements :
- GATK 4.2.5.0
- Samtools 1.15
- bwa 0.7.17
- Python 3.10.2
- Cromwell 76 (using WDL 1.0)
- Docker support

### Input requirements/expectations:
- Cleaned-up paired-end reads in FASTQ format and/or unmapped BAM or CRAM files for `initial` mode (one pair/file per sample).
- Mapped and de-duplicated BAM or CRAM files, with indexes, for `repeat` and `final` modes.
- Read group information, either provided with the JSON input for FASTQ files, or included in the uBAM or BAM files. 
- A reference genome indexed using either bwa or bwa-mem2 (support for bwa-mem2 is planned for the longer term).
- Truth SNP and/or Indel sets, if available, in either uncompressed or gzip-compressed VCF format with TBI index(es).
- A list of 'scatters', either produced with ScatterIntervalsByNs (to be described in forthcoming documentation) or with the third-party Python tool, [ScaffoldStitcher](https://github.com/ameliahaj/ScaffoldStitcher). This facilitates scatter-gathering across the reference genome.
- If in `repeat` or `final` mode, the `polymorphicRegions.tar.gz` file produced in `initial` mode.

## Configuring the JSON input file:

The pipeline is optimized for user inputs to be as simple and limited as possible. Beyond the required files (_i.e._ reference and index files; FASTQ, BAM or uBAM files; outputs from prior modes), all configuration and parameters are supplied in a single JSON input file.

### Mandatory options

**The following options must be set in the input JSON file:**

| Mode | Option | Type | Description |
| --- | ------ | ---- | ---- |
| All | `mode` | String | `initial`, `repeat`, or `final` |
| All | `ref` | File | Reference file in `.vcf` or `.vcf.gz` format |
| All | `ref_dict` | File | Reference dictionary file in `.dict` format |
| All | `ref_fai` | File | Reference index file in `.fai` format |
| Initial | `bwamem2` | Boolean | Default is false. Set 'true' to use bwa-mem2 instead of bwa mem. |
| Initial | `ref_amb` | File | Required if using bwa mem or bwa-mem2 |
| Initial | `ref_ann` | File | Required if using bwa mem or bwa-mem2 |
| Initial | `ref_pac` | File | Required if using bwa mem or bwa-mem2 |
| Initial | `ref_bwt` | File | Required if using bwa mem or bwa-mem2 |
| Initial | `ref_sa` | File | Required if using bwa mem only |
| Initial | `ref_0123` | File | Required if using bwa-mem2 only |
| Initial | `ref_bwt_2bit_64` | File | Required if using bwa-mem2 only |
| Repeat/Final | `packaged_polymorphic_regions` | File | Output from the prior run in `initial` mode. |
<br />

**The following must always be provided for each user-defined scatter, as a `scatterList` object:**

| Name | Type | Description |
| ------ | ---- | ---- |
| `name` | String | Name of user-defined scatter |
| `intervals` | String | Comma-separated list in chr:pos format, _e.g._ `chr21:1-27800000` |
<br />

**The following must always be provided for each sample, as a `sampleList` object:**

| Name | Type | Description |
| ------ | ---- | ---- |
| `name` | String | Name of sample or individual. |
| `taxon_group` | String | Taxonomic group name. This is used to group individuals by taxon. |
<br />

**The `sampleList` object must also contain the following, depending on mode and sample input type:**

| Input Type | Name | Mode | Type |  Description |
| ------ | ------ | ---- | ---- | ---- |
| FASTQ: | `R1` | Initial | File | First of the paired-end FASTQ files.  |
| - | `R2` | Initial | File | Second of the paired-end FASTQ files.  |
| - | `RG_ID` | Initial | String | Read group ID.  |
| - | `RG_SM` | Initial | String | Read group sample name. |
| - | `RG_PU` | Initial | String | Read group platform unit. Note this is used in BQSR, which models together all reads with the same PU. |
| uBAM: | `unmapped_bam` | Initial | File | Unmapped BAM file if not mapping from FASTQ; must contain all read group data.  |
| BAM: | `bam` | Repeat/Final | File | Recalibrated BAM file, _e.g._ from previous mode, required in each of these modes. |
| - | `bam_index` | Repeat/Final | File |  Recalibrated BAM file index, required in each of these modes.  |

### Common optional inputs

**The following options may also be set, depending on mode of operation:**

| Mode | Option | Type | Description |
| --- | --- | --- | --- |
| Initial | `cram_not_bam` | Boolean | Default is true; indicates that CRAM versus BAM format should be used throughout the workflow. |
| Initial | `flowcell_patterned` | Boolean | Default is true; sets `--optical-duplicate-pixel-distance` in `MarkDuplicatesSpark` to 2500. Set `false` for 100. |
| Any | `truth_set_SNPs` | File? | Truth set of known SNPs in `.vcf` or `.vcf.gz` format. |
| Any | `truth_set_SNPs_index` | File? | Index to the above file, in `.tbi` format. |
| Any | `truth_set_INDELs` | File? | Truth set of known Indels in `.vcf` or `.vcf.gz` format. |
| Any | `truth_set_INDELs_index` | File? | Index to the above file, in `.tbi` format. |
| Any | `validate_truth_sets` | Boolean | Set `true` to perform `ValidateVariants` on the truth sets. |
| Any | `merge_contigs_into_num_partitions` | Int? | Optional parameter for `GenomicsDBImport`. |
| Any | `container_gatk` | String? | Default is `broadinstitute/gatk:4.2.0.0`. |
| Any | `container_gitc` | String? | Default is `us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.5.7-2021-06-09_16-47-48Z`. |
| Any | `container_python` | String? | Default is `python:3.9.5`. |
| Any | `path_to_gitc` | String? | Default is `/usr/gitc`. |
| Any | `path_to_gitc_gatk` | String? | Default is `/usr/gitc/gatk4/`. |
<br />
Other options are available for advanced users; specifically, to configure the resources assigned to container instances. If you know how to use these correctly, you'll also know how to find them!

## Important notes and caveats to consider

- The pipeline was originally written in the WDL Development Spec for use with Cromwell 65. However, because Terra supports only WDL 1.0 at present, the current version has been re-written to be 'backwards compatible' with WDL 1.0. For this reason, the code may not flow intuitively in places (_e.g._ some WDL dev functions are written as 'pseudo-functions' in separate WDLs for 1.0: see `/workflows/functions/CollectByKey.wdl` for an example). The upside is that these can be simply switched out when WDL is upgraded.
- Allele-specific annotations are requested throughout, though full support of [allele-specific filtering](https://gatk.broadinstitute.org/hc/en-us/articles/360035890551?id=9622) is not yet implemented.
- If no truth sets are provided, or if only one truth set is provided, multi-allelic variants are first split to be biallelic (using `LeftAlignAndTrimVariants`) and then split into separate SNP and Indel files (using `SplitVcfs`). This is imperative, else — at least in allele-specific mode — [mixed sites will fail to be processed correctly](https://gatk.broadinstitute.org/hc/en-us/articles/360035890551?id=9622). If both truth sets are provided, the unfiltered variants are recalibrated with `ApplyVQSR` run in series, as recommended.
- The --merge-contigs-into-num-partitions parameter is set to 0 (default), but can (and probably should) be configured to a higher value in the input JSON [as explained here](https://gatk.broadinstitute.org/hc/en-us/articles/360056138571-GDBI-usage-and-performance-guidelines). We are still determining best practices for handling this when large numbers of contigs are in use, or when scatters/interval lists split scaffolds into smaller contigs. Note that this parameter requires whole contigs versus intervals.
- Re-genotyping at the final callset capitalizes on the `--include-non-variant-sites` option in `GenotypeGVCFs` [as implemented in GATK 4.0.12.0](https://github.com/broadinstitute/gatk/releases/tag/4.0.12.0).
- The pipeline currently uses the Broad's production container for 'Genomes in the Cloud', which uses outdated versions of all tools (_i.e._ GATK 4.1.8.0 versus 4.2.5.0; bwa 0.7.15 versus 0.7.17; samtools 1.11 versus 1.15). This container is only used when mapping reads or when running `AnalyzeCovariates`. Eventually, we'll replace this with a new container comprising the latest required software versions ([per this issue](https://github.com/broadinstitute/GATK-For-Primates/issues/9)).
- When running `AnalyzeCovariates`, the required 'R' libraries are downloaded and installed each time. This is not optimal!

## Acknowledgements

This pipeline is a collaborative effort by the [Broad Institute](https://www.broadinstitute.org) ([@bhanugandham](https://github.com/bhanugandham)) and ([@grahamlbanes](https://github.com/grahamlbanes)), with support from the [Institute of Museum and Library Services](https://www.imls.gov) and [Google Cloud for Education](https://cloud.google.com/edu/). Development was partially facilitated by the compute resources and assistance of the [Center For High Throughput Computing (CHTC)](https://chtc.cs.wisc.edu) in the Department of Computer Sciences at the [University of Wisconsin-Madison](https://www.wisc.edu/). The CHTC is supported by the [Advanced Computing Initiative](https://aci.wisc.edu), the [Wisconsin Alumni Research Foundation](https://www.warf.org/), the [Wisconsin Institutes for Discovery](https://discovery.wisc.edu), and the [National Science Foundation (NSF)](https://www.nsf.gov/), and is an active member of the [Open Science Grid](https://opensciencegrid.org), which is supported by NSF and the [U.S. Department of Energy's Office of Science](https://www.energy.gov/science/office-science).
