#!/bin/bash -ue
bwa mem -R "@RG\tID:Seq01p\tSM:Seq01\tPL:ILLUMINA\tPI:330" -t 10 -M human_g1k_v37_decoy.small.fasta input.1 C09DF_N_111207.4.AGTTGCTT_R1_xxx.fastq.gz > 1234N.sam
