#!/bin/bash

declare -a pipeline_steps_somatic=("make_bam_and_qc.nf", "somatic.nf")

for i in "$@"
do
case $i in
    -n=*|--nextflow=*)
    NEXTFLOW="${i#*=}"
    shift # past argument=value
    ;;
    -a=*|--analysis=*)
    ANALYSIS="${i#*=}"
    shift # past argument=value
    ;;
    -m=*|--mapping=*)
    MAPPING="${i#*=}"
    shift # past argument=value
    ;;
    -p=*|--pairing=*)
    PAIRING="${i#*=}"
    shift # past argument=value
    ;;
    -pr=*|--profile=*)
    PROFILE="${i#*=}"
    shift # past argument=value
    ;;
    *)
    echo "Unknown option $i" # unknown option
    exit 1
    ;;
esac

done
if [[ -n $1 ]]; then
    echo "Last line of file specified as non-opt/last argument:"
    tail -1 $1
fi

if [ $ANALYSIS = "somatic" ]; 
then
    eval $NEXTFLOW run make_bam_and_qc.nf --mapping $MAPPING --pairing $PAIRING -profile $PROFILE && $NEXTFLOW run somatic.nf --sample make_bam_output.tsv -profile $PROFILE
fi
