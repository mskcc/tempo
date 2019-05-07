#!/bin/bash

NEXTFLOW=""
ANALYSIS=""
MAPPING=""
PAIRING=""
PROFILE=""


function help () {
  echo "Usage:"
  echo "    run_pipeline.sh --nextflow=</path/to/nextflow> --analysis=<analysis> --mapping=/path/to/mapping_file.tsv --pairing=/path/to/pairing_file.tsv"
  echo "Options:"
  echo "    -n --nextflow: path to nextflow executable"
  echo "    -a --analysis: somatic or germline"
  echo "    -m --mapping: mapping file"
  echo "    -p --pairing: pairing file"
  echo "    -pr --profile: profile"
  exit 1
}

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
    help
    ;;
esac
done



if [ "$NEXTFLOW" = "" ] || [ "$ANALYSIS" = "" ] || [ "$MAPPING" = "" ] || [ "$PAIRING" = "" ] || [ "$PROFILE" = "" ];
then
    help
fi

if [[ -n $1 ]]; then
    echo "Last line of file specified as non-opt/last argument:"
    tail -1 $1
fi

if [ "$ANALYSIS" = "make_bam" ]; 
then
    eval $NEXTFLOW run make_bam_and_qc.nf --mapping $MAPPING --pairing $PAIRING -profile $PROFILE
fi

if [ "$ANALYSIS" = "somatic" ]; 
then
    eval $NEXTFLOW run make_bam_and_qc.nf --mapping $MAPPING --pairing $PAIRING -profile $PROFILE && $NEXTFLOW run somatic.nf --sample make_bam_output.tsv -profile $PROFILE
fi

if [ "$ANALYSIS" = "germline" ];
then
    eval $NEXTFLOW run make_bam_and_qc.nf --mapping $MAPPING --pairing $PAIRING -profile $PROFILE && $NEXTFLOW run germline.nf --sample make_bam_output.tsv -profile $PROFILE
fi
