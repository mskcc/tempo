#!/bin/bash

nextflow run main_align_markDups_BaseRecal.nf --sample ./test_samples.tsv --outDir s3://testbucketcmo/nextflow-test -profile awsbatch
