# Running the Pipeline

## Overview

::: tip Note
Tempo does not support running samples from mixed sequencing platforms together. By default, the pipeline assumes the inputs are from exome sequencing.
:::

This page provides instructions on how to run the pipeline through the `pipeline.nf` script. The basic command below shows how to run Tempo, with an explanation of flags and input arguments and files. Below is also described how to best [run the pipeline on Juno](running-the-pipeline.md#running-the-pipeline-on-juno) as well as [on AWS](running-the-pipeline.md#running-the-pipeline-on-aws).

```shell
nextflow run pipeline.nf \
    --mapping/--bamMapping <input mapping tsv file> \
    --pairing <input pairing tsv file, can be optional> \
    --assayType <string value: either "exome" or "genome"> \
    --outDir <path to output subdirectory> \ 
    -profile juno \
    --somatic --germline --QC\
    --aggregate <true, false, or [a tsv file]>
```

_Note: [The number of dashes matters](nextflow-basics.md)._

**Required arguments:**
* `--mapping/--bamMapping <tsv>` is required except running in `--aggregate [tsv]` mode. When `--mapping [tsv]` is provided, FASTQ file paths are expected in the TSV file, and the pipeline will start from FASTQ files and go through all steps to generate BAM files. When `--bamMapping [tsv]` if provided, BAM file paths are expected in the TSV file. See [The Mapping File](running-the-pipeline.md#input-files) and [Execution Mode](running-the-pipeline.md#execution-mode) for details.
* `--pairing <tsv>` is required when `--somatic` and/or `--germline` are enabled. `--pairing <tsv>` is not needed when you are running BAM generation part only, even if you are doign it with `--QC`( or `--QC` and `--aggregate`) enabled. See [The Mapping File](running-the-pipeline.md#input-files) and [Execution Mode](running-the-pipeline.md#execution-mode) for details.
* `--assayType` ensures appropriate resources are allocated for indicated assay type. Only `exome` or `genome` is supported. Note: Please also make sure this value matches the `TARGET` field you put in the mapping.tsv file. Available TARGET field value for `exome` are `idt` or `agilent`i (can be mixed), for `genome` is `wgs`.
* `--outDir` is the directory where the output will end up. This directory does not need to exist. If not set, by default it will be set to run directory (i.e. the directory from which the command `nextflow run` is executed.)
* `-profile` loads the preset configuration required to run the pipeline in the supported environment. Accepted values are `juno` and `awsbatch` for execution on the [Juno cluster](juno-setup.md) or on [AWS Batch](aws-setup.md), respectively. `-profile test_singularity` is for testing on `juno`.

**Section arguments:**
* `--somatic`, `--germline` and `--QC` flags are boolean that indicate to run the somatic, germline variant calling and QC modules, respectively. Default value are `false` for all. Note: Currently the pipeline will enable `--somatic` automatically if only `--germline` is enabled, since germline analysis needs results from somatic analysis for now. (default: `false` for all)
* `--aggregate <true, false, or [a tsv file]>` can be boolean or be given a path to a tsv file. Default value is `false`. It has to work together with `--somatic`, `--germline` and/or `--QC` to aggregate the results of these operations together as a cohort. There will be a `cohort_level/[cohort]/` directory generated under `--outDir [path]`. When boolean value `true` is given (equal to only give `--aggregate`), TEMPO will aggregate all the samples in the mapping and pairing file as one cohort named "default cohort". When `--aggregate <tsv>` file is given, the pipeline will aggregate samples and tumor/normal pairs based on the value is given in `COHORT` columns. Each sample and tumor/normal pairs can be assigned to different cohorts in different rows.


**Optional arguments:**
* `-work-dir`/`-w` is the directory where the temporary output will be cached. By default, this is set to the run directory. Please see `NXF_WORK` in [Nextflow environment variables](https://www.nextflow.io/docs/latest/config.html#environment-variables).
* `-publishAll` is a boolean, resulting in retention of intermediate output files ((default: `true`).
* `--splitLanes` indicates that the provided FASTQ files will be scanned for all unique sequencing lanes and demultiplexed accordingly. This is recommended for some steps of the alignment pipeline. See more under [The Mapping File](running-the-pipeline.md#input-files) (default: `true`).
* `-with-timeline` and `-with-report` are enabled by default and results in the generation of a timeline and resource usage report for the pipeline run. These are boolean but can also be fed output names for the respective file.
* `--genome` is the version of reference files for your analysis. Currently only `GRCh37` is supported. We will add support for `GRCh38` later. (default: `GRCh37`)

Using test inputs provided in the GitHub repository, here is a concrete example:

```shell
nextflow run pipeline.nf \
    --mapping test_inputs/local/full_test_mapping.tsv \ 
    --pairing test_inputs/local/full_test_pairing.tsv \
    -profile juno \
    --outDir results
    --somatic --germline --aggregate --QC \
```

## Input Files

::: tip Note
The header lines are mandatory in the following files, but not the order of their columns. The nesessary header fields must be included, and additional columns are allowed.
:::

::: warning Be aware
Tempo checks for the following aspects:
1. Duplicated rows
2. Valid headers
3. Concordance of given `--assayType` and `TARGET` field in mapping file. (for `exome` are `idt` or `agilent`, for `genome` is `wgs`)
4. Concordance of `TARGET` field in mapping file for Tumor/Normal pairs that is defined in pairing file.
5. File paths in the mapping file are all valid
6. File extentions in the mapping file
7. BAI files exisits in the same directory for the BAM files in the mapping file. (Basename can be either `*.bai` or `*.bam.bai`)
8. Required Arguments are provided and valid
9. Reference files all exists
10. Samples are present both in mapping and pairing files at the same time
11. Duplicated files are in the mapping file
But we suggest you do your own validation of your input to ensure smooth execution.
:::

### The Mapping Files

#### FASTQ Mapping File (`--mapping <tsv>`)

For processing paired-end FASTQ inputs, users must provide both a mapping file using `--mapping <tsv>`, as described below.

You do not need `--pairing <tsv>` when you are running BAM generation alone (even together with `--QC`).

You must to give `--pariring <tsv>` when `--somatic` or `--germline` is enabled.

::: warning Be aware
Tempo can deal with any number of sequencing lanes per sample, in any combination of lanes split or combined across multiple FASTQ pairs. Different FASTQ pairs for the same sample can be provided as different lines in the mapping file and give the same SAMPLE ID in the `SAMPLE` field, and repeating the `TARGET` field. By default, Tempo will look for all distinct sequencing lanes in provided FASTQ files by scanning each FASTQ read name. The pipeline uses this and the instrument, run, and flowcell IDs from the _sequence identifiers_ in the input FASTQs to generate all different read group IDs for each sample. This information is used by the base quality score recalibration steps of the GATK suite of tools. If FASTQ files name explicitly specified the lane name in the format of `_L(\d){3}_` ("\_L" + "3 integer" + "\_"), the pipeline will assume this FASTQ files contain only one lane, and it will skip scanning and splitting the FASTQ files, and give one read group ID all the reads in the FASTQ files based on the name of the first read in the FASTQ file. Please refer to [this GATK Forum Page](https://gatkforums.broadinstitute.org/gatk/discussion/6472/read-groups)for more details.
:::

This file is necessary to map the input FASTQ pairs from one or more FASTQ pairs to SAMPLE IDs. Additionally, this file tells the pipeline what bait set it used (`wgs` in the case of whole genome sequencing).

Example:

|SAMPLE|TARGET|FASTQ_PE1|FASTQ_PE2|
|:---:|:---:|:---:|:---:|:---:|
|normal_sample_1|agilent|normal1_L001_R01.fastq.gz|normal1_L001_R02.fastq.gz|
|normal_sample_1|agilent|normal1_L002_R01.fastq.gz|normal1_L002_R02.fastq.gz|
|tumor_sample_1|agilent|tumor1_L001_R01.fastq.gz|tumor1_L001_R02.fastq.gz|
|tumor_sample_1|agilent|...|...|
|tumor_sample_1|agilent|tumor1_L00N_R01.fastq.gz|tumor1_L00N_R02.fastq.gz|

Accepted values for the **TARGET** column are `agilent`, `idt` or `wgs`. Please note `idt` and `agilent` can be mixed and are valid when `--assayType exome`. But `wgs` can not be mixed with any other value, and is only valid when `--assayType genome`\

Read further details on these parameters [here](reference-files.md#genomic-intervals).

#### BAM Mapping File (`--bamMapping <tsv>`)

If the user is using pre-processed BAMs, the input TSV file is a similar format as FASTQ mapping TSV file, with slight difference showing below.

You must to give `--pariring <tsv>` and at least one of `--somatic` `--germline` `--QC` simultaneously when you use `--bamMapping <tsv>`.

Example:

|SAMPLE|TARGET|BAM|BAI|
|:---:|:---:|:---:|:---:|:---:|
|normal_sample_1|agilent|normal1.bam|normal1.bai|
|normal_sample_2|agilent|normal2.bam|normal2.bai|
|tumor_sample_1|agilent|...|...|
|tumor_sample_2|agilent|tumor2.bam|tumor2.bai|

The `--pairing <tsv>` file will be exactly the same as using FASTQ mapping TSV file, describing below.

::: tip Note
The pipeline expects BAM file indices in the same subdirectories as `TUMOR_BAM` and `NORMAL_BAM`. If the index files `*.bai` or `*.bam.bai` do not exist, `pipeline.nf` will throw an error. The BAI column in the BAM Mapping TSV file is not actually used.

Different from FASTQ mapping tsv, in this TSV file each SAMPLE id can only appear once, meaning the pipeline will not combine different BAMs for you for the same sample.
:::

### The Pairing File (`--pairing <tsv>`)

The pipeline needs to know which tumor and normal samples are to be analyzed as matched pairs. This file provides that pairing by referring to the sample names as provided in the `SAMPLE` column in the mapping file.

You do not need `--pairing <tsv>` when you are running BAM generation alone (even together with `--QC`).

Example:

|NORMAL_ID|TUMOR_ID|
|:---:|:---:|
|normal_sample_1|tumor_sample_1|
|normal_sample_2|tumor_sample_2|
|...|...|
|normal_sample_n|tumor_sample_n|

### Aggregate File (`--aggregate true/false/<tsv>`)

* When boolean value `true` is given (equal to only give `--aggregate`), TEMPO will aggregate all the samples in the mapping and pairing file as one cohort named "default cohort". 
* When `--aggregate <tsv>` file is given, the pipeline will aggregate samples and tumor/normal pairs based on the value is given in `COHORT` column. Each sample and tumor/normal pairs can be assigned to different cohorts in different rows.
* When running aggregation only mode, `PATH` column need to be provided to introduce the TEMPO result directories for each sample and tumor/normal pairs (only up to the parent folder of `qc`, `somatic` and `germline` folder).

Example:

|NORMAL_ID|TUMOR_ID|COHORT|PATH(only applicable for runing aggregate only mode)|
|:---:|:---:|:---:|:---:|:---:|
|normal_sample_1|normal_sample_1|cohort1|/home/tempo/v1/result|
|normal_sample_2|tumor_sample_2|cohort1|/home/tempo/v1/result|
|normal_sample_3|tumor_sample_3|...|...|
|normal_sample_n|tumor_sample_n|cohort2|/home/tempo/v2/result|

The `--pairing <tsv>` file will be exactly the same as using FASTQ mapping TSV file, describing below.

::: tip Note
The pipeline expects BAM file indices in the same subdirectories as `TUMOR_BAM` and `NORMAL_BAM`. If the index files `*.bai` or `*.bam.bai` do not exist, `pipeline.nf` will throw an error. The BAI column in the BAM Mapping TSV file is not actually used.

Different from FASTQ mapping tsv, in this TSV file each SAMPLE id can only appear once, meaning the pipeline will not combine different BAMs for you for the same sample.
:::

## Execution Mode

There are overall 4 execution modes that TEMPO accept, depending on the way of input arguments are given. The pipeline will throw an error when it detects incompatible _Section arguments_ (`--somatic`, `--germline`, `--QC`, `--aggregate`) combinations. Compatibility of analysis argument combinations are described below:

### `--mapping <tsv>` only
When no additional _Section arguments_ are given, only alignment steps will be performed.

***Compatible _Section Arguments_ Combinations:***
* `--QC`: `QcAlfred` and `QcCollectHsMetrics` will be performed.
* `--QC --aggregate true/<tsv>`: `QcAlfred` and `QcCollectHsMetrics` will be performed and aggregated in `cohort_level/[cohort]` folder.

***Incompatible _Section Arguments_:***
* `--somatic`: No pairing infomation.
* `--germline`: No pairing information.
* `--bamMapping <tsv>`: Conflicts.

### `--mapping/--bamMapping <tsv>` and `--pairing <tsv>` (We are describing two modes in this section)
When `--mapping <tsv>` is given, the pipeline will use FASTQ input and start from alignment steps.

When `--bamMapping <tsv>` is given, the pipeline will use BAM input and skip alignment steps.

When no additional _Sectionarguments_ are given, pipeline will throw an error indicating that `--pairing <tsv>` is not used.

***Compatible _Section Arguments_ Combinations:***
* `--somatic` with or without `--aggregate true/<tsv>`
* `--somatic --germline` with or without `--aggregate true/<tsv>`
* `--somatic --QC` with or without `--aggregate true/<tsv>`
* `--QC` with or without `--aggregate true/<tsv>`: `QcConpair` will be performed together with `QcAlfred` and `QcCollectHsMetrics`.

***Incompatible _Section Arguments_:***
* `--germline` with or without `--aggregate true/<tsv>` and `--QC`: The pipeline will auto-enable `--somatic` because germline analysis need the results from somatic analysis at this stage.
* `--mapping <tsv>`: Conflicts.

### `--aggregate <tsv>` only
This mode can only be run when the TEMPO produced output structure path is provided as `PATH` column in `--aggregate <tsv>`. It explicitly relies on the output structure ((only up to the parent folder of `qc`, `somatic` and `germline` folder) that are auto-generated by TEMPO to identify how and what files need to be aggregated together as a cohort level result under folder `cohort_level/[cohort]`. Please refer to [Outputs](outputs.md#outputs) for more detail.


***Compatible _Section Arguments_ Combinations:***
* The pipeline will auto-detect which sections need to be aggregated, so no _Section_Arguments_ need to be given.

***Incompatible _Section Arguments_:***
* `--somatic`: Not needed. Auto-detected.
* `--germline`: Not needed. Auto-detected.
* `--QC`: Not needed. Auto-detected.
* `--mapping <tsv>`: Conflicts.
* `--bamMapping <tsv>`: Conflicts.
* `--pairing <tsv>`: Conflicts.

## Running the Pipeline on Juno

::: tip Note
First follow the instructions to [set up your enviroment on Juno](juno-setup.md).
:::

### Submitting the Pipeline to LSF

We recommend submitting your `nextflow run pipeline.nf <...>` command to the cluster via `bsub`, which will launch a leader job from which individual processes are submitted as jobs to the cluster.

```
bsub -W <hh:mm> -n 2 -R "rusage[mem=<requested memory>]" \
    -o <LSF output file name>.out -e <LSF error file name>.err \
    nextflow run pipeline.nf -profile juno <...> 
```

We recommend that users check the [documentation for LSF](https://www.ibm.com/support/knowledgecenter/en/SSETD4_9.1.2/lsf_command_ref/bsub.1.html) to clarify each of the arguments above. However,

* `-W <hh:mm>` sets the time allotted for `nextflow run pipeline.nf` to run to completion. 
* `-n 2` is requesting one slot. This should be sufficient for `nextflow run pipeline.nf`
* ` -o <LSF output file name>.out` is the name of the STDOUT file, which is quite informative for Nextflow. We **strongly** encourage users to set this.
* ` -e <LSF output file name>.err` is the name of the STDERR file. Please set this. 
* ` -R "rusage[mem=<requested memory>]"` is the requested memory for  `nextflow run pipeline.nf`, which will not be memory intensive at all. 

Here is a concrete example of a bsub command to process 25 WES TN pairs, running somatic and germline variant calling modules:

```shell
bsub -W 80:00 -n 2 -R "rusage[mem=8]" -o nf_output.out -e nf_output.err \
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

It is normally not a good idea to run things on the log-in nodes of the cluster. Instead we recommend scheduling an interactive session via e.g. ` bsub -Is -n 2 -R "rusage[mem=8]" csh` and running the `screen` within that session.

Users are welcome to use `nohup` or `tmux` as well. 

## Running the Pipeline on AWS

::: tip Note
These instructions will assume the user is moderately knowledgeable of AWS. Please refer to [AWS Setup](aws-setup.md) and the [AWS Glossary](aws-glossary.md) we have curated.
:::

## Modifying or Resuming Pipeline Run

Nextflow supports [modify and resume](https://www.nextflow.io/docs/latest/getstarted.html?#modify-and-resume). 

To resume an interrupted Nextflow pipeline run, add `-resume` (note the single dash) to your command-line call to access the cache history of Nextflow and continue a job from where it left off. This will trigger a check of which jobs already completed before starting unfinished jobs in the pipeline.

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
