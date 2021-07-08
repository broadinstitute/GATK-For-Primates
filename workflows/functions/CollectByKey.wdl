version 1.0

## Copyright Broad Institute and Wisconsin National Primate Research Center,
## University of Wisconsin-Madison, 2021
## 
## Functions in lieu of WDL 1.1/dev SPEC functions that are not yet
## implemented in Terra/Cromwell, from the complete germline short
## variant discovery pipeline optimized for non-human primates.
## For requirements, expectations and outputs, please review the
## complete documentation at:
## https://github.com/broadinstitute/GATK-For-Primates
##
## LICENSING : 
## This script is released under the WDL source code license (BSD-3) (see LICENSE in 
## https://github.com/broadinstitute/wdl). However, the programs it calls may be
## subject to different licenses. Users are responsible for checking that they are
## authorized to run all programs before running this script. Please see the docker
## containers for detailed licensing information pertaining to the included programs.

##########################################################################
## WORKFLOW DEFINITION
##########################################################################

workflow collectByKey {

    input {
        Array[String] allGroups
        Array[String] groups
        Array[String] members
    }

        scatter (y in allGroups) {
            scatter (x in zip(groups,members)) {
                if (x.left == y) {
                    String? z = x.right
                    String? z_index = x.right + ".tbi"
                }
                
            }
            Pair[String,Array[String]] main = (y, select_all(z))
            Pair[String,Array[String]] main_index = (y, select_all(z_index))
        }

    output {
        Array[Pair[String,Array[String]]] collected = main
        Array[Pair[String,Array[String]]] collected_indexes = main_index
    }

}
