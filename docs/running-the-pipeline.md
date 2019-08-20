# Running the Pipeline

## Overview 

This page provides instructions on how to run the pipeline through the `pipeline.nf` script. The basic command below shows how to run Vaporware, with an explanation of flags and input arguments and files. Below is also described how to best [run the pipeline on Juno](run-pipeline.md#running-the-pipeline-on-juno) as well as [on AWS](run-pipeline.md#running-the-pipeline-on-aws).

```shell
nextflow run pipeline.nf \
    --somatic --germline --outDir <path to output subdirectory> \ 
    -w <path to temporary work directory with cached results> \
    -profile juno \
    --mapping <input mapping tsv file>  \
    --pairing <input pairing tsv file>  \
    -assayType <string value: either "exome" or "genome">
    -with-report <name of html file> \ 
    -with-trace <name of html file> 
```
**NOTE:** The number of dashes [matters](nextflow-basics.md)

* The `--somatic` and `--germline` flags indicate to run the somatic and germline variant calling modules, respectively. If not set, `pipeline.nf` will only align BAMs.
* `--outDir` is the directory where the output will end up. This directory does not need to exist. If not set, by default it will be set to run directory (i.e. the directory from which the command `nextflow run` is executed.)
* `-w` is the directory where the temporary output will be cached. By default, this is set to the run directory. Please see `NXF_WORK` in [Nextflow environment variables](https://www.nextflow.io/docs/latest/config.html#environment-variables).
* `-profile` loads the preset configuration required to run the pipeline in the supported environment. Accepted values are `juno` and `awsbatch` for execution on the [Juno cluster](juno-setup.md) or on [AWS Batch](aws-setup.md), respectively.
* The files provided to the `--mapping` and `--pairing` arguments should contain the mapping of FASTQ files to sample names and of tumor-normal pairs. These are tab-separated files, see further description below and examples in the `test_inputs` subdirectory.

Using test inputs provided in the GitHub repository, here is a concrete example:

```shell
nextflow run pipeline.nf --somatic --germline \
    --mapping test_inputs/local/full_test_mapping.tsv \ 
    --pairing test_inputs/local/full_test_pairing.tsv \
    -profile juno \
    --outDir results
```

## Input Files

For processing paired-end FASTQ inputs, users must provide both a mapping file and pairing file, as described below.

### The Mapping File

This file is necessary to map the input FASTQ pairs from one or more sequencing lanes to sample names. Additionally, this file tells the pipeline whether the samples are exome or genome samples. In the case of the former, the capture kit used is also input.

_Note: The header line is mandatory but not the order of columns._

Example:

|SAMPLE|LANE|ASSAY|TARGET|FASTQ_PE1|FASTQ_PE2|
|:---:|:---:|:---:|:---:|:---:|:---:|
|normal_sample_1|L001|wes|agilent|normal1_L001_R01.fastq.gz|normal1_L001_R02.fastq.gz|
|normal_sample_1|L002|wes|agilent|normal1_L002_R01.fastq.gz|normal1_L002_R02.fastq.gz|
|tumor_sample_1|L001|wes|agilent|tumor1_L001_R01.fastq.gz|tumor1_L001_R02.fastq.gz|
|tumor_sample_1|...|wes|agilent|...|...|
|tumor_sample_1|L00N|wes|agilent|tumor1_L00N_R01.fastq.gz|tumor1_L00N_R02.fastq.gz|

Accepted values for the **ASSAY** column are `exome` and `genome`.\
Accepted values for the **TARGET** column are `agilent` and `idt`.\
Read further details on these parameters [here](bioinformatics-components.md#genome-versus-exome).

_Note: It is important that **LANE** values be descriptive, as the combination of *SAMPLE* and *LANE* values must be unique. Duplicate non-unique values cause errors._

### The Pairing File
The pipeline needs to know which tumor and normal samples are to be analyzed as matched pairs. This files provides that pairing by referring to the sample names as provided in the **SAMPLE** column in the mapping file.

_Note: The header line is mandatory but not the order of columns._

Example:

|NORMAL_ID|TUMOR_ID|
|:---:|:---:|
|normal_sample_1|tumor_sample_1|
|normal_sample_2|tumor_sample_2|
|...|...|
|normal_sample_n|tumor_sample_n|

## Running with Input BAMs

If the user is processing input BAMs, a mapping of tumor and normal sample names to BAM files and their pairing is used as below, in lieu of providing separate mapping and pairing files as when starting from FASTQ files. 

```shell
nextflow run pipeline.nf --somatic --germline \
    --bam_pairing <input bam TN pairs tsv file> \
    -profile juno \
    --outDir results 
```
The `bam_pairing` input file is also a tab-separated file, see a further description below.

### The BAM Pairing File
Given BAMs as inputs, the user must specify which tumor and normal samples are to be analyzed as matched pairs. The following format is used:

_Note 1: The header line is mandatory but not the order of columns._
_NOTE 2: The pipeline expects BAM file indices in the same subdirectories as `TUMOR_BAM` and `NORMAL_BAM`. If the index files `*.bai` do not exist, `pipeline.nf` will throw an error requiring users to do so._

Example:

|TUMOR_ID|NORMAL_ID|ASSAY|TARGET|TUMOR_BAM|NORMAL_BAM|
|:---:|:---:|:---:|:---:|:---:|:---:|
|normal_sample_1|tumor_sample_1|wes|agilent|/path/to/file/tumor_1.bam|/path/to/file/normal_1.bam|
|normal_sample_2|tumor_sample_2|wes|agilent|/path/to/file/tumor_2.bam|/path/to/file/normal_2.bam|
|...|...|...|...|...|...|
|normal_sample_n|tumor_sample_n|wes|agilent|/path/to/file/tumor_n.bam|/path/to/file/normal_n.bam|

## Running the Pipeline on Juno

### Submitting the Pipeline to LSF

We recommend submitting your `nextflow run pipeline.nf <...>` command to the cluster via `bsub`, which will launch a leader job from which individual processes are submitted as jobs to the cluster.

```
bsub -W <hh:mm> -n 1 -R "rusage[mem=<requested memory>]" \
    -o <LSF output file name>.out -e <LSF error file name>.err \
    nextflow run pipeline.nf -profile juno <...> 
```

It is **important** that users use a recent version of singularity, as detailed in [Juno setup](juno-setup.md).

We recommend that users check the [documentation for LSF](https://www.ibm.com/support/knowledgecenter/en/SSETD4_9.1.2/lsf_command_ref/bsub.1.html) to clarify each of the arguments above. However,

* `-W <hh:mm>` sets the time allotted for `nextflow run pipeline.nf` to run to completion. 
* `-n 1` is requesting one slot. This should be sufficient for `nextflow run pipeline.nf`
* ` -o <LSF output file name>.out` is the name of the STDOUT file, which is quite informative for Nextflow. We **strongly** encourage users to set this.
* ` -e <LSF output file name>.err` is the name of the STDERR file. Please set this. 
* ` -R "rusage[mem=<requested memory>]"` is the requested memory for  `nextflow run pipeline.nf`, which will not be memory intensive at all. 


Here is a concrete example of a bsub command to process 25 WES TN pairs, running somatic and germline variant calling modules:

```shell
bsub -W 80:00 -n 1 -R "rusage[mem=15]" -o nf_output.out -e nf_output.err \
    nextflow run <path-to-repository>/pipeline.nf --somatic --germline \
    --mapping test_inputs/local/WES_25TN.tsv --pairing test_inputs/local/WES_25TN_pairing.tsv 
    --outDir results \
    -profile juno
```

### Running From a `screen` Session

Another option is to use a `screen` session for running the pipeline interactively, for example naming and entering a screen session as follows:

`screen -RD new_screen_name`

It is normally not a good idea to run things on the log-in nodes of the cluster. Instead we recommend scheduling an interactive session via e.g. ` bsub -Is -n 1 -R "rusage[mem=20]" csh` and running the `screen` within that session.

Users are welcome to use `nohup` or `tmux` as well. 

## Running the Pipeline on AWS

These instructions will assume the user is moderately knowledgeable of AWS. Please refer to [AWS Setup](aws-setup.md) and the [AWS Glossary](aws-glossary.md) we have curated.

## Modifying or Resuming Pipeline Run

Nextflow supports [modify and resume](https://www.nextflow.io/docs/latest/getstarted.html?#modify-and-resume). 

To resume an interrupted Nextflow pipeline run, add `-resume` (note the single dash) to your command-line call to access Nextflow's cache history and continue a job from where it left off. This will trigger a check of which jobs already completed before starting unfinished jobs in the pipeline.

This function also allows you to make changes to values in the `pipeline.nf` script and continue from where you left off. Nextflow will use the cached information from the unchanged sections while running only the modified processes. If you want to make changes to processes that already successfully completed, you have to manually delete the subdirectories in `work` where those processes where run. 

_Note 1: if you use `-resume` for the first time of a timeline run, Nextflow will recognize this as superfluous, and continue._

_Note 2: To peacefully interrupt an ongoing Nextflow pipeline run, do `control+C` once and wait for Nextflow to kill submitted jobs. Otherwise orphan jobs might be left on the cluster._