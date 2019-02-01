#!/bin/bash -ue
gatk ApplyBQSR   -R human_g1k_v37_decoy.small.fasta   --input 1234N_XX.md.bam   --output 1234N.recal.bam   -L small.intervals   --create-output-bam-index true   --bqsr-recal-file 1234N.recal.table
