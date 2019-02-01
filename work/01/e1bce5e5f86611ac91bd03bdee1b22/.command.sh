#!/bin/bash -ue
samtools sort -m 2G -@ 10 -m 2G -o 1234N.sorted.bam 1234N.bam
