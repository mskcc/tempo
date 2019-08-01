# Run the Vaporware pipeline

This basic command shows how to run Vaporware, using test input provided in the GitHub repository (in the `REPO` directory). Below the input arguments and files are explained:
```shell
nextflow run pipeline.nf \
    --somatic --germline --outDir test_output -profile juno \
    --mapping $REPO/test_inputs/local/test_make_bam_and_qc.tsv \
    --pairing $REPO/test_inputs/local/test_make_bam_and_qc_pairing.tsv
```

* The `--somatic` and `--germline` flags indicate to run the somatic and germline variant calling modules, respectively. 
* `--outDir` is the directory where the output will end up. This directory does not need to exist. If not set, it will be the directory from which the command is run.
* `-profile` loads the preset configuration required to run the pipeline in the supported environment. Accepted values are `juno` and `awsbatch` for execution on the [Juno cluster](juno-setup.md) or on [AWS Batch](aws-setup.md).
* The files provided to the `--mapping` and `--pairing` arguments should contain the mapping of FASTQ files to sample names and of tumor-normal pairs. These are tab-separated files, see further description below and examples in the [test_inputs](test_inputs) directory.

#### The mapping file
This file is necessary to map the input FASTQ pairs from one or more sequencing lanes to sample names. Additionally, this file tells the pipeline whether the samples are exome or genome samples. In the case of the former, the capture kit used is also input.

_Note: The header line is mandatory but not the order of columns._

|SAMPLE|LANE|ASSAY|TARGET|FASTQ_PE1|FASTQ_PE2|
|:---:|:---:|:---:|:---:|:---:|:---:|
|normal_sample_1|L001|wes|agilent|L001_R01.fastq.gz|L001_R02.fastq.gz|
|normal_sample_1|L002|wes|agilent|L002_R01.fastq.gz|L002_R02.fastq.gz|
|tumor_sample_1|L001|wes|agilent|L001_R01.fastq.gz|L001_R02.fastq.gz|
|tumor_sample_1|...|wes|agilent|...|...|
|tumor_sample_1|L00N|wes|agilent|L00N_R01.fastq.gz|L00N_R02.fastq.gz|

Accepted values for the **ASSAY** column are `exome` and `genome`.\
Accepted values for the **TARGET** column are `agilent` and `idt`.\
Read further details on these parameters [here](bioinformatics-components.md#genome-versus-exome).

#### The pairing file
The pipeline needs to know which tumor and normal samples are to be analyzed as matched pairs. This files provides that pairing by referring to the sample names as provided in the **SAMPLE** column in the mapping file.

_Note: The header line is mandatory but not the order of columns._

Example:
|NORMAL_ID|TUMOR_ID|
|:---:|:---:|
|normal_sample_1|tumor_sample_1|
|normal_sample_2|tumor_sample_2|
|...|...|
|normal_sample_n|tumor_sample_n|

### Running with pre-aligned BAM files