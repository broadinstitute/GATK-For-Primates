version development

## Adapted from the generate-sample-map.wdl (C) Broad Institute 2020
## https://github.com/gatk-workflows/utility-wdls/blob/main/generate-sample-map.wdl

task BAMs {

    input {
        Array[String] bams
        Array[String] bam_groups
        Array[String] bam_names
        Array[String] bam_indexes
        String docker = "python:latest"
    }

    command <<<
    set -oe pipefail
    
    python << CODE
    bams = ['~{sep="','" bams}']
    bam_groups = ['~{sep="','" bam_groups}']    
    bam_names = ['~{sep="','" bam_names}']
    bam_indexes = ['~{sep="','" bam_indexes}']

    if len(bams)!= len(bam_indexes) != len(bam_names) != len(bam_groups):
      print("Numbers of input variables are not equal.")
      exit(1)

    with open("bams_generated_tsv.txt", "w") as fi:
      for i in range(len(bams)):
        fi.write(bam_names[i] + "\t" + bam_groups[i] + "\t" + bams[i] + "\t" + bam_indexes[i] + "\n")
    CODE

    tsv_of_bams = read_tsv(bams_generated_tsv.txt)
    >>>

    runtime {
        docker: docker
    }

    output {
        File tsv = "bams_generated_tsv.txt"
        Array[Array[String]] tsv_of_bams = tsv_of_bams
    }

}