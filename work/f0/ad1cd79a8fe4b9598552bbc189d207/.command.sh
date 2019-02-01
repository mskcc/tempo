#!/bin/bash -ue
gatk MarkDuplicates --java-options "-Xms4000m -Xmx7g"    --MAX_RECORDS_IN_RAM 50000   --INPUT 1234N.sorted.bam   --METRICS_FILE 1234N.bam.metrics   --TMP_DIR .   --ASSUME_SORT_ORDER coordinate   --CREATE_INDEX true   --OUTPUT 1234N_XX.md.bam
