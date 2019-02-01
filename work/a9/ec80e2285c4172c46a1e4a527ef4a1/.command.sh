#!/bin/bash -ue
gatk MarkDuplicates --java-options "-Xms4000m -Xmx7g"    --MAX_RECORDS_IN_RAM 50000   --INPUT 9876T.sorted.bam   --METRICS_FILE 9876T.bam.metrics   --TMP_DIR .   --ASSUME_SORT_ORDER coordinate   --CREATE_INDEX true   --OUTPUT 9876T_XX.md.bam
