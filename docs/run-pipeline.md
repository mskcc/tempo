# Running the Vaporware pipeline

## Overview 

This basic command shows how to run Vaporware, with an explanation of flags and input arguments:
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
* `--outDir` is the directory where the output will end up. This directory does not need to exist. If not set, by default it will be set to run directory (i.e. the directory from which 'the command `nextflow run` is executed.)
* `-w` is the directory where the temporary output will be cached. By default, this is set to the run directory. Please see `NXF_WORK` in [Nextflow environment variables](https://www.nextflow.io/docs/latest/config.html#environment-variables)
* `-profile` loads the preset configuration required to run the pipeline in the supported environment. Accepted values are `juno` and `awsbatch` for execution on the [Juno cluster](juno-setup.md) or on [AWS Batch](aws-setup.md).
* The files provided to the `--mapping` and `--pairing` arguments should contain the mapping of FASTQ files to sample names and of tumor-normal pairs. These are tab-separated files, see further description below and examples in the test_inputs subdirectory.
* `-assayType` loads the standard configured resources to run either WES or WGS via Juno or AWS. Accepted values are the string `"exome"` (which is default) or `"genome"`. If the flag is left used, it's assumed the user is processing exomes.

Using test input provided in the GitHub repository (in the `REPO` directory), here is a concrete example:

```
nextflow run pipeline.nf --somatic --germline \
    --mapping test_inputs/local/full_test_mapping.tsv \ 
    --pairing test_inputs/local/full_test_pairing.tsv \
     -profile juno --outDir Result \
     -with-report report.html
     -with-timeline timeline.html
```

### Running with Input BAMs

If the user is processing input BAMs, use the following command:
```shell
nextflow run pipeline.nf \
    --somatic --germline --outDir <path to output subdirectory> \ 
    -w <path to temporary work directory with cached results> \
    -profile juno \
    --bam_pairing <input bam TN pairs tsv file>  \
    -assayType <string value: either "exome" or "genome">
    -with-report <name of html file> \ 
    -with-trace <name of html file> 
```

Note the only difference with the above is that users are required to input a TN pairing file of BAMs via the argument `--bam_pairing`. Do not use the arguments `--mapping` or `--pairing` in this case, as this is reserved for processing FASTQs. The `bam_pairing` input file is also a tab-separated `*txt` or `*tsv` file; please see a further description below and examples in the test_inputs subdirectory in the Vaporware repo. 

## Input files

For processing paired-end FASTQ inputs, users must provide both a mapping file and pairing file:

### The mapping file
This file is necessary to map the input FASTQ pairs from one or more sequencing lanes to sample names. Additionally, this file tells the pipeline whether the samples are exome or genome samples. In the case of the former, the capture kit used is also input.

_Note: The header line is mandatory but not the order of columns._

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

**NOTE** It is very important that **LANE** values be descriptive, as the combination of SAMPLE_LANE values must be unique. Duplicate non-unique values cause errors.

For instance, 

|SAMPLE|LANE|ASSAY|TARGET|FASTQ_PE1|FASTQ_PE2|
|:---:|:---:|:---:|:---:|:---:|:---:|
|normal_sample_1|8|wes|agilent|normal_sample_1_AA.R1.fastq.gz|normal_sample_1_AA.R2.fastq.gz|
|normal_sample_1|8|wes|agilent|normal_sample_1_BB.R1.fastq.gz|normal_sample_1_BB.R2.fastq.gz|


will cause a fatal error, as there are two rows with `normal_sample_1_8`, even if FASTQ filenames differ. 


### The pairing file
The pipeline needs to know which tumor and normal samples are to be analyzed as matched pairs. This files provides that pairing by referring to the sample names as provided in the **SAMPLE** column in the mapping file.

_Note: The header line is mandatory but not the order of columns._

Example:

|NORMAL_ID|TUMOR_ID|
|:---:|:---:|
|normal_sample_1|tumor_sample_1|
|normal_sample_2|tumor_sample_2|
|...|...|
|normal_sample_n|tumor_sample_n|



### The BAM pairing file (for BAM inputs)
Given BAMs as inputs, the user must specify which tumor and normal samples are to be analyzed as matched pairs. The following format is used:

_Note: The header line is mandatory but not the order of columns._

Example:

|TUMOR_ID|NORMAL_ID|ASSAY|TARGET|TUMOR_BAM|NORMAL_BAM|
|:---:|:---:|:---:|:---:|:---:|:---:|
|normal_sample_1|tumor_sample_1|wes|agilent|/path/to/file/tumor_1.bam|/path/to/file/normal_1.bam|
|normal_sample_2|tumor_sample_2|wes|agilent|/path/to/file/tumor_2.bam|/path/to/file/normal_2.bam|
|...|...|...|...|...|...|
|normal_sample_n|tumor_sample_n|wes|agilent|/path/to/file/tumor_n.bam|/path/to/file/normal_n.bam|


**NOTE:** We expect that users have indexed BAMs with `samtools index` in the same subdirectories as the BAMs in `TUMOR_BAM` and `NORMAL_BAM`. If the BAM index files `*.bai` do not exist, `pipeline.nf` will throw an error requiring users to do this step. 



## Running the pipeline on Juno

### LSF

We recommend submitting `nextflow run pipelin.nf` via bsub:

```
bsub -W <hh:mm> -n 1 -R "rusage[mem=<requested memory>]" \
    -o <LSF output file name>.out -e <LSF error file name>.err nextflow run pipeline.nf --somatic --germline \
    --mapping <input mapping tsv file> --pairing <input pairing tsv file> -profile juno \
    --outDir <path to output subdirectory> -w <path to temporary work directory with cached results> \
    -with-report <name of html report file> -with-trace <name of html trace file> 
```

It's **important** that users use the recent version of singularity, as detailed in [Juno Setup](juno-setup.md).

Before submitting this bsub command, either load the relevant module via 

```shell
module load singularity/3.1.1
```

or export the correct version of singularity via 

```shell
export PATH=/opt/local/singularity/3.1.1/bin:$PATH
```

We recommend that users check the [documentation for LSF](https://www.ibm.com/support/knowledgecenter/en/SSETD4_9.1.2/lsf_command_ref/bsub.1.html) to clarify each of the arguments above. However,

* `-W <hh:mm>` sets the time alloted for `nextflow run pipeline.nf` to run to completion. 
* `-n 1` is requesting one slot. This should be sufficient for `nextflow run pipeline.nf`
* ` -o <LSF output file name>.out` is the name of the STDOUT file, which is quite informative for Nextflow. We **strongly** encourage users to set this.
* ` -e <LSF output file name>.err` is the name of the STDERR file. Please set this. 
* ` -R "rusage[mem=<requested memory>]"` is the requested memory for  `nextflow run pipeline.nf`, which will not be memory intensive at all. 


Here is a concrete example of a bsub command to process 25 WES TN pairs, both somatic and germline:

```
bsub -W 80:00 -n 1 -R "rusage[mem=15]" -o nf_output.out -e nf_output.err /juno/work/taylorlab/biederstedte/vaporware/nextflow run /juno/work/taylorlab/biederstedte/vaporware/pipeline.nf -w /juno/work/taylorlab/biederstedte/vaporware/work  --mapping /juno/work/taylorlab/biederstedte/vaporware/test_inputs/local/WES_25TN.tsv  --pairing /juno/work/taylorlab/biederstedte/vaporware/test_inputs/local/WES_25TN_pairing.tsv  --somatic --germline --outDir /juno/work/taylorlab/biederstedte/vaporware/output -profile juno  -with-report report_25WESTN.html -with-trace -with-timeline timeline_25WESTN.html
```


### screen or tmux

Another option is to use a screen session for running the pipeline interactively, naming and entering a screen session as follows:

`$ screen -RD new_screen_name`

Once you have set up the appropriate variables for singularity and the Nextflow singularity cache with `export NXF_SINGULARITY_CACHEDIR` (see the details in [Juno setup](juno-setup.md)), run the pipeline as follows: 

```
nextflow run pipeline.nf \
    --somatic --germline --outDir <path_to_output_subdirectory> -w <path_to_temp_work_subdirectory> -profile juno \
    --mapping <path_to_pairing_file>
    --pairing <path_to_mapping_file>
```

*NOTE:* Normally it's not a good idea to run things on the log-in nodes of the server. Consider scheduling an intervative sessioon via e.g. ` bsub -Is -n 1 -R "rusage[mem=20]" csh`

Users welcome to use **nohup** or **tmux** as well. 


## Running the pipeline on AWS

These instructions will assume the user is moderately knowlecgable of AWS. Please refer to [AWS Setup](aws-setup.md) and the [AWS Glossary](aws-glossary.md) we have curated.






