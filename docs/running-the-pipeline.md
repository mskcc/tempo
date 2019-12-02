# Running the Pipeline

## Overview 

::: tip Note
Tempo does not support running samples from mixed sequencing platforms together.
:::

This page provides instructions on how to run the pipeline through the `pipeline.nf` script. The basic command below shows how to run Tempo, with an explanation of flags and input arguments and files. Below is also described how to best [run the pipeline on Juno](running-the-pipeline.md#running-the-pipeline-on-juno) as well as [on AWS](running-the-pipeline.md#running-the-pipeline-on-aws).

```shell
nextflow run pipeline.nf \
    --somatic --germline \
    --assayType <string value: either "exome" or "genome"> \
    --outDir <path to output subdirectory> \ 
    -profile juno \
    --mapping <input mapping tsv file> \
    --pairing <input pairing tsv file>
```

_Note: [The number of dashes matters](nextflow-basics.md)._

**Recommended arguments:**
* The `--somatic` and `--germline` flags are boolean that indicate to run the somatic and germline variant calling modules, respectively. If not set, the pipeline will only align BAMs.
* `--assayType` ensures appropriate resources are allocated for indicated assay type.
* `--outDir` is the directory where the output will end up. This directory does not need to exist. If not set, by default it will be set to run directory (i.e. the directory from which the command `nextflow run` is executed.)
* `-profile` loads the preset configuration required to run the pipeline in the supported environment. Accepted values are `juno` and `awsbatch` for execution on the [Juno cluster](juno-setup.md) or on [AWS Batch](aws-setup.md), respectively.
* The files provided to the `--mapping` and `--pairing` arguments should contain the mapping of FASTQ files to sample names and of tumor-normal pairs. These are tab-separated files, see further description below and examples in the [test inputs subdirectory](https://github.com/mskcc/tempo/tree/master/test_inputs).

**Optional arguments:**
* `-work-dir`/`-w` is the directory where the temporary output will be cached. By default, this is set to the run directory. Please see `NXF_WORK` in [Nextflow environment variables](https://www.nextflow.io/docs/latest/config.html#environment-variables).
* `-publishAll` is a boolean, resulting in retention of intermediate output files ((default: `true`).
* `--splitLanes` indicates that the provided FASTQ files will be scanned for all unique sequencing lanes and demultiplexed accordingly. This is recommended for some steps of the alignment pipeline. See more under [The Mapping File](running-the-pipeline.md#input-files) (default: `true`).
* `-with-timeline` and `-with-report` are enabled by default and results in the generation of a timeline and resource usage report for the pipeline run. These are boolean but can also be fed output names for the respective file.
* `--conpairAll` runs the Conpair sample concordance assessment for all combinations of tumor and normal samples in the run (default: `false`).

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

::: tip Note
The header lines are mandatory in the following files, but not the order of their columns.
:::

::: warning Be aware
Tempo checks for duplicated combinations of sample and lane names, empty entries, and some other things. However, it is up to the user to make sure that the inputs are in good shape.
:::

::: warning Be aware
Tempo can deal with any number of sequencing lanes per sample, in any combination of lanes split or combined across multiple FASTQ pairs. By default, Tempo will look for all distinct sequencing lanes in provided FASTQ files. The pipeline uses this and the instrument, run, and flowcell IDs from the _sequence identifiers_ in the input FASTQs to generate all different read group IDs for each sample. This information is used by the base quality score recalibration steps of the GATK suite of tools. If FASTQ files name explicitly specified the lane name in the format of `_L(\d){3}_` ("_L" + "3 integer" + "_"), the pipeline will assume this FASTQ files contain only one lane, and it will skip scanning and splitting the FASTQ files, and give one read group ID all the reads in the FASTQ files.
:::

### The Mapping File

This file is necessary to map the input FASTQ pairs from one or more sequencing lanes to sample names. Additionally, this file tells the pipeline whether the samples are exome or genome samples. In the case of the former, the capture kit used is also input.

Example:

|SAMPLE|ASSAY|TARGET|FASTQ_PE1|FASTQ_PE2|
|:---:|:---:|:---:|:---:|:---:|
|normal_sample_1|wes|agilent|normal1_L001_R01.fastq.gz|normal1_L001_R02.fastq.gz|
|normal_sample_1|wes|agilent|normal1_L002_R01.fastq.gz|normal1_L002_R02.fastq.gz|
|tumor_sample_1|wes|agilent|tumor1_L001_R01.fastq.gz|tumor1_L001_R02.fastq.gz|
|tumor_sample_1|wes|agilent|...|...|
|tumor_sample_1|wes|agilent|tumor1_L00N_R01.fastq.gz|tumor1_L00N_R02.fastq.gz|

Accepted values for the **ASSAY** column are `exome` and `genome`.\
Accepted values for the **TARGET** column are `agilent` and `idt`.\
Read further details on these parameters [here](reference-resources.md#genomic-intervals).

### The Pairing File

The pipeline needs to know which tumor and normal samples are to be analyzed as matched pairs. This files provides that pairing by referring to the sample names as provided in the **SAMPLE** column in the mapping file.

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
    --bamPairing <input bam TN pairs tsv file> \
    -profile juno \
    --outDir results 
```
The `bamPairing` input file is also a tab-separated file, see a further description below.

### The BAM Pairing File
Given BAMs as inputs, the user must specify which tumor and normal samples are to be analyzed as matched pairs. The following format is used:

Example:

|TUMOR_ID|NORMAL_ID|ASSAY|TARGET|TUMOR_BAM|NORMAL_BAM|
|:---:|:---:|:---:|:---:|:---:|:---:|
|normal_sample_1|tumor_sample_1|wes|agilent|/path/to/file/tumor_1.bam|/path/to/file/normal_1.bam|
|normal_sample_2|tumor_sample_2|wes|agilent|/path/to/file/tumor_2.bam|/path/to/file/normal_2.bam|
|...|...|...|...|...|...|
|normal_sample_n|tumor_sample_n|wes|agilent|/path/to/file/tumor_n.bam|/path/to/file/normal_n.bam|

::: tip Note
The pipeline expects BAM file indices in the same subdirectories as `TUMOR_BAM` and `NORMAL_BAM`. If the index files `*.bai` do not exist, `pipeline.nf` will throw an error.
:::

## Running the Pipeline on Juno

::: tip Note
First follow the instructions to [set up your enviroment on Juno](juno-setup.md).
:::

### Submitting the Pipeline to LSF

We recommend submitting your `nextflow run pipeline.nf <...>` command to the cluster via `bsub`, which will launch a leader job from which individual processes are submitted as jobs to the cluster.

```
bsub -W <hh:mm> -n 1 -R "rusage[mem=<requested memory>]" \
    -o <LSF output file name>.out -e <LSF error file name>.err \
    nextflow run pipeline.nf -profile juno <...> 
```

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

::: warning Be aware
Whereas a few exome samples finish within a few hours, larger batches and genomes will take .s. Allow for this by setting`-W` to a good amount of hours. The pipeline will die if the leader job does, but can be [resumed](running-the-pipelinf.md#modifying-or-resuming-pipeline-run) subsequently. 
:::

### Running From a `screen` Session

Another option is to use a `screen` session for running the pipeline interactively, for example naming and entering a screen session as follows:

`screen -RD new_screen_name`

It is normally not a good idea to run things on the log-in nodes of the cluster. Instead we recommend scheduling an interactive session via e.g. ` bsub -Is -n 1 -R "rusage[mem=20]" csh` and running the `screen` within that session.

Users are welcome to use `nohup` or `tmux` as well. 

## Running the Pipeline on AWS

::: tip Note
These instructions will assume the user is moderately knowledgeable of AWS. Please refer to [AWS Setup](aws-setup.md) and the [AWS Glossary](aws-glossary.md) we have curated.
:::

*Under development.*

## Modifying or Resuming Pipeline Run

Nextflow supports [modify and resume](https://www.nextflow.io/docs/latest/getstarted.html?#modify-and-resume). 

To resume an interrupted Nextflow pipeline run, add `-resume` (note the single dash) to your command-line call to access Nextflow's cache history and continue a job from where it left off. This will trigger a check of which jobs already completed before starting unfinished jobs in the pipeline.

This function also allows you to make changes to values in the `pipeline.nf` script and continue from where you left off. Nextflow will use the cached information from the unchanged sections while running only the modified processes. If you want to make changes to processes that already successfully completed, you have to manually delete the subdirectories in `work` where those processes where run. 

::: tip Note
* If you use `-resume` for the first time of a timeline run, Nextflow will recognize this as superfluous, and continue.
* To peacefully interrupt an ongoing Nextflow pipeline run, do `control+C` once and wait for Nextflow to kill submitted jobs. Otherwise orphan jobs might be left on the cluster.
:::

To resume the pipeline from a specific run, please read the pages here on using [resume](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html)and as well troubleshooting [resumed runs](https://www.nextflow.io/blog/2019/troubleshooting-nextflow-resume.html) for more complicated use cases.

In order to resume from a specific time you ran the pipeline, first check the specific pipeline runs with `nextflow log`:

```shell
> nextflow log

TIMESTAMP            DURATION  RUN NAME          STATUS  REVISION ID  SESSION ID                            COMMAND                                    
2019-05-06 12:07:32  1.2s      focused_carson    ERR     a9012339ce   7363b3f0-09ac-495b-a947-28cf430d0b85  nextflow run hello                         
2019-05-06 12:08:33  21.1s     mighty_boyd       OK      a9012339ce   7363b3f0-09ac-495b-a947-28cf430d0b85  nextflow run rnaseq-nf -with-docker        
2019-05-06 12:31:15  1.2s      insane_celsius    ERR     b9aefc67b4   4dc656d2-c410-44c8-bc32-7dd0ea87bebf  nextflow run rnaseq-nf                     
2019-05-06 12:31:24  17s       stupefied_euclid  OK      b9aefc67b4   4dc656d2-c410-44c8-bc32-7dd0ea87bebf  nextflow run rnaseq-nf -resume -with-docker
```

Users can then restart the pipeline at specific run, using either the `RUN NAME` or the `SESSION ID`. For instance

```shell
> nextflow run rnaseq-nf -resume mighty_boyd
```

or equivalently

```shell
> nextflow run naseq-nf -resume 4dc656d2-c410-44c8-bc32-7dd0ea87bebf
```

Sometimes the resume feature may not work entirely as expected, as described in troubleshooting tips [here on the Nextflow blog](https://www.nextflow.io/blog/2019/troubleshooting-nextflow-resume.html)


## After Successful Run

Nextflow generate many intermediate output files. All the relevant output data should be in the directory given to the `outDir` argument. Once you have verified that the data are satisfactory, everything outside this directory can be removed. In particular, the `work` directory will contain all intermediate output files, which takes up a great deal of disk space and should be removed. The `nextflow clean -force` command does all of this. Also see `nextflow clean -help` for options. 

::: warning Be aware
Once these files are removed, modifications to or resumption of a pipeline run **cannot** be done.
:::
