# Vaporware 

CMO WES/WGS pipeline using [Nextflow framework](https://github.com/nextflow-io/nextflow)

See Nextflow documention here: 
https://www.nextflow.io/

Inspiration from `Sarek` at the NBIS and SciLifeLab in Sweden:
https://github.com/SciLifeLab/Sarek

Vaporware contains several pipeline `*.nf` scripts but all functionality is contained in `pipeline.nf`.

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

## Description of pipeline.nf

The `pipeline.nf` script takes a mapping and pairing file to perform BAM preprocessing, somatic analysis, and germline analysis. BAM Preprocessing always occurs; to activate the somatic and germline pipelines, use the `--somatic` and `--germline`, respectively, at the command line. See the sections below for command examples.

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

The default parameters are for local-use *WITHOUT* containers.

* For Docker use, do the following:

```
nextflow run pipeline.nf --somatic --germline --mapping mapping.tsv --pairing pairing.tsv -profile docker
```

* For Singularity use, do the following:

```
nextflow run pipeline.nf --somatic --germline --mapping mapping.tsv --pairing pairing.tsv -profile singularity
```

## Useful Tips

* __Number of dashes matter__: Nextflow has a quirk where Nextflow executor-specific, built-in flags are initiated with a single dash, whereas parameters we define require two. For example, to define which run profile to use, we call Nextflow's built-in feature `-profile`; to do resume, we use `-resume`. To provide mapping or pairing input file paths, which call functions the development team at MSKCC had to write, we have to call it with `--mapping` and `--pairing`, respectively.

* __Modify and Resume__: Nextflow supports [modify and resume](https://www.nextflow.io/docs/latest/getstarted.html?#modify-and-resume); to activate this feature, use `-resume` (note the single dash) at the command line to use Nextflow's cache history and continue a job from where it left off. This even lets you make changes to values in the script and continue from there - Nextflow will use the cached information from the unchanged sections while running only the modified processes.

* __View the Nextflow Log__: You can access the Nextflow cache metadata by running `nextflow log`. It contains information like TIMESTAMP, DURATION, and RUN_NAME, and a STATUS indicating a failed or successful run, among others. The values under RUN_NAME can also be submitted following the `-resume` flag to resume previously-run Nextflow jobs.

* __Run/Skip specific tools:__ `pipeline.nf` lets you run only specific tools you want and skip others. The aforementioned `--somatic` and `--germline` flags already have a preset of tools to run, but you can limit this further by providing the `--tools` flag, followed by a comma-delimited string. For example, to use only DELLY for your somatic/germline runs, do `--somatic --germline --tools delly`; to use MuTect2, Manta, and Strelka2, do `--somatic --tools mutect2,manta,strelka2`.

## Components 

The following links to the wiki provide greater details:

* [Bioinformatic components](https://github.com/mskcc/vaporware/wiki/Bioinformatic-Components)
* [Building an AWS Batch Environment](https://github.com/mskcc/vaporware/wiki/Building-AWS-Batch-Compute-Environment)
* [Reference Files](https://github.com/mskcc/vaporware/wiki/Reference-Files)
* [Genomic Intervals](https://github.com/mskcc/vaporware/wiki/Genomic-Intervals)
* [Variant Annotation](https://github.com/mskcc/vaporware/wiki/Variant-Annotation)


