#!/bin/bash -ue
gatk BaseRecalibrator   --input 9876T_XX.md.bam   --output 9876T.recal.table   --tmp-dir /tmp   -R human_g1k_v37_decoy.small.fasta   -L small.intervals   --known-sites dbsnp_138.b37.small.vcf   --known-sites Mills_and_1000G_gold_standard.indels.b37.small.vcf --known-sites 1000G_phase1.indels.b37.small.vcf   --verbosity INFO
