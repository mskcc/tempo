#!/bin/bash

trap _term SIGINT

NEXTFLOW=""
ANALYSIS=""
MAPPING=""
PAIRING=""
PROFILE=""


_term() {
  echo "Killing nextflow process" 
  kill -TERM "$make_bam_pid" 2>/dev/null
  kill -TERM "$somatic_pid" 2>/dev/null
  kill -TERM "$germline_pid" 2>/dev/null
}


function help () {
  echo "Usage:"
  echo "    run_pipeline.sh --nextflow=</path/to/nextflow> --analysis=<analysis> --mapping=/path/to/mapping_file.tsv --pairing=/path/to/pairing_file.tsv"
  echo "Options:"
  echo "    -n --nextflow: path to nextflow executable"
  echo "    -a --analysis: make_bam, somatic, germline"
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

FAIL=0

if [ "$NEXTFLOW" = "" ] || [ "$ANALYSIS" = "" ] || [ "$MAPPING" = "" ] || [ "$PAIRING" = "" ] || [ "$PROFILE" = "" ];
then
    help
fi

if [[ -n $1 ]]; then
    echo "Last line of file specified as non-opt/last argument:"
    tail -1 $1
fi

echo "RUNNING with" $($NEXTFLOW -v)

if [ "$ANALYSIS" = "make_bam" ]; 
then
    echo "Starting MAKE_BAM_AND_QC analysis"
    echo "Starting make_bam_and_qc.nf step"
    eval $NEXTFLOW run make_bam_and_qc.nf --mapping $MAPPING --pairing $PAIRING -profile $PROFILE &
    make_bam_pid=$!
    echo "Started make_bam_and_qc.nf. You can track progress in make_bam_and_qc.log"
    wait $make_bam_pid || let "FAIL+=1"
    if [ "$FAIL" == "0" ];
    then
        echo "Step make_bam_and_qc COMPLETED"
    else
        echo "STEP: make_bam_and_qc FAILED. You can see log in make_bam_and_qc.log"
        exit 1
    fi
fi

if [ "$ANALYSIS" = "germline" ]; 
then
    echo "Starting GERMLINE analysis"
    echo "Starting make_bam_and_qc.nf step"
    $NEXTFLOW run make_bam_and_qc.nf --mapping $MAPPING --pairing $PAIRING -profile $PROFILE > make_bam_and_qc.log &
    make_bam_pid=$!
    echo "Started make_bam_and_qc.nf. You can track progress in make_bam_and_qc.log"
    wait $make_bam_pid || let "FAIL+=1"
    if [ "$FAIL" == "0" ];
    then
        echo "Step make_bam_and_qc COMPLETED"
        echo "Starting somatic.nf step"
        $NEXTFLOW run somatic.nf --sample make_bam_output.tsv -profile $PROFILE > somatic.log &
        somatic_pid=$!
        echo "Started somatic.nf. You can track progress in somatic.log"
        echo "Starting germline.nf step"
        $NEXTFLOW run germline.nf --sample make_bam_output.tsv -profile $PROFILE > germline.log &
        germline_pid=$!
        echo "Started germline.nf. You can track progress in germline.log"
        wait $somatic_pid || let "FAIL+=1"
        if [ "$FAIL" == "0" ];
        then
            echo "Step somatic COMPLETED"
        else
            echo "Step somatic FAILED. You can see log in somatic.log"
            exit 1
        fi
        wait $germline_pid || let "FAIL+=1"
        if [ "$FAIL" == "0" ];
        then
            echo "Step germline COMPLETED"
        else
            echo "Step germline FAILED. You can see log in make_bam_and_qc.log"
            exit 1
        fi
    else
        echo "STEP: make_bam_and_qc FAILED. You can see log in make_bam_and_qc.log"
        exit 1
    fi
fi

if [ "$ANALYSIS" = "somatic" ];
then
    echo "Starting SOMATIC analysis"
    echo "Starting make_bam_and_qc.nf step"
    $NEXTFLOW run make_bam_and_qc.nf --mapping $MAPPING --pairing $PAIRING -profile $PROFILE > make_bam_and_qc.log &
    make_bam_pid=$!
    echo "Started make_bam_and_qc.nf. You can track progress in make_bam_and_qc.log"
    wait $make_bam_pid || let "FAIL+=1"
    if [ "$FAIL" == "0" ];
    then
        echo "Step make_bam_and_qc COMPLETED"
        echo "Starting somatic.nf step"
        $NEXTFLOW run somatic.nf --sample make_bam_output.tsv -profile $PROFILE > somatic.log &
        somatic_pid=$!
        echo "Started somatic.nf. You can track progress in somatic.log"
        wait $somatic_pid || let "FAIL+=1"
        if [ "$FAIL" == "0" ];
        then
            echo "Step somatic COMPLETED"
        else
            echo "Step somatic FAILED. You can see log in somatic.log"
            exit 1
        fi
    else
        echo "STEP: make_bam_and_qc FAILED. You can see log in make_bam_and_qc.log"
        exit 1
    fi
fi
