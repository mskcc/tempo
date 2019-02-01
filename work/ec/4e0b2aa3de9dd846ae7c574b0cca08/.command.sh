#!/bin/bash -ue
gatk ApplyBQSR   -R human_g1k_v37_decoy.small.fasta   --input 9876T_XX.md.bam   --output 9876T.recal.bam   -L small.intervals   --create-output-bam-index true   --bqsr-recal-file 9876T.recal.table
