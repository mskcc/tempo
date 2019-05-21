# Vaporware 

CMO WES/WGS pipeline using [Nextflow framework](https://github.com/nextflow-io/nextflow)

See Nextflow documention here: 
https://www.nextflow.io/

Inspiration from `Sarek` at the NBIS and SciLifeLab in Sweden:
https://github.com/SciLifeLab/Sarek

Vaporware contains several pipeline scripts (`pipeline.nf`, `make_bam_and_qc.nf`, `somatic.nf`, `germline.nf`). A detalied explanation for each step is provided below.

## Setup instructions

Clone the repository as `vaporware`:

```
git clone https://github.com/mskcc/vaporware.git vaporware
cd vaporware
```

Install `Nextflow` within this subdirectory if the Nextflow executable `nextflow` isn't in your PATH:

```
curl -s https://get.nextflow.io | bash 
```

This installs Nextflow within the subdirectory `vaporware`.

Optionally, move the `nextflow` file to a directory accessible by your `$PATH` variable in order to avoid remembering and typing the full path to `nextflow` each time.


## For submitting via LSF on juno

We recommend users use `singularity/3.1.1`, i.e. 

```
$ module load  singularity/3.1.1
$ which singularity
/opt/local/singularity/3.1.1/bin/singularity
```

**NOTE:** In order to run successfully on `juno`, you must set `NXF_SINGULARITY_CACHEDIR`. Use e.g. 
`export NXF_SINGULARITY_CACHEDIR=$HOME/.singularity`

See here for details: https://www.nextflow.io/docs/latest/config.html

```
NXF_SINGULARITY_CACHEDIR:
Directory where remote Singularity images are stored. 
When using a computing cluster it must be a shared folder accessible from all computing nodes.
```

* Do the following for LSF on juno:

```
nextflow run pipeline.nf --somatic --germline --mapping test_inputs/lsf/test_make_bam_and_qc.tsv --pairing test_inputs/lsf/test_make_bam_and_qc_pairing.tsv -profile juno
```

## For submitting via AWS Batch

In order to run pipeline on `AWS Batch`, you first must create your `Compute Environment` and `Configuration File` as described [here](aws_cf_scripts/README.md).
 
* When you build your compute environment and create configuration file do the following:

```
nextflow run pipeline.nf --somatic --germline --mapping test_inputs/aws/test_make_bam_and_qc.tsv --pairing test_inputs/aws/test_make_bam_and_qc_pairing.tsv -profile awsbatch
```

## Docker and Singularity

**NOTE:** To being able to run locally you need to provide reference files from `conf/references.config` and create `samples.tsv` as described in the wiki page [Bioinformatic components](https://github.com/mskcc/vaporware/wiki/Bioinformatic-Components)

The default parameters are for local-use *WITHOUT* containers

* For Docker use, do the following:

```
nextflow run pipeline.nf --somatic --germline --mapping mapping.tsv --pairing pairing.tsv -profile docker
```

* For Singularity use, do the following:

```
nextflow run pipeline.nf --somatic --germline --mapping mapping.tsv --pairing pairing.tsv -profile singularity
```

## Components 

The following links to the wiki provide greater details:

* [Bioinformatic components](https://github.com/mskcc/vaporware/wiki/Bioinformatic-Components)
* [Building an AWS Batch Environment](https://github.com/mskcc/vaporware/wiki/Building-AWS-Batch-Compute-Environment)
* [Reference Files](https://github.com/mskcc/vaporware/wiki/Reference-Files)
* [Genomic Intervals](https://github.com/mskcc/vaporware/wiki/Genomic-Intervals)
* [Variant Annotation](https://github.com/mskcc/vaporware/wiki/Variant-Annotation)


