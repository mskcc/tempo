#!/bin/bash -ue
samtools sort -m 2G -@ 10 -m 2G -o 9876T.sorted.bam 9876T.bam
