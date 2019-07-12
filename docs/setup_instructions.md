## Setup Instructions

### Using vaporware

Clone the repository:

```
git clone https://github.com/mskcc/vaporware.git 
cd vaporware
```

**NOTE** Users can name this repo anything via `git clone https://github.com/mskcc/vaporware.git differentNameHere`


Next, users must install `Nextflow` within this subdirectory if the Nextflow executable `nextflow` isn't in your PATH:

```
curl -s https://get.nextflow.io | bash 
```

This installs Nextflow within the subdirectory `vaporware`.

Optionally, move the `nextflow` file to a directory accessible by your `$PATH` variable in order to avoid remembering and typing the full path to `nextflow` each time.


### Vaporeware configuration scripts

Nextflow allows users to set up 

See here for details: https://www.nextflow.io/docs/latest/config.html



## For submitting via LSF on juno

We recommend users use `singularity/3.1.1`, i.e. 

```
$ module load singularity/3.1.1
$ which singularity
/opt/local/singularity/3.1.1/bin/singularity
```

Users are recommended to set the following environmental variable to call this version of singularity by default: 

```
export PATH=/opt/local/singularity/3.1.1/bin:$PATH
```

**NOTE:** In order to run successfully on `juno`, you must set `NXF_SINGULARITY_CACHEDIR`. Use e.g. 
`export NXF_SINGULARITY_CACHEDIR=$HOME/.singularity`

See here for details: https://www.nextflow.io/docs/latest/config.html

```
NXF_SINGULARITY_CACHEDIR:
Directory where remote Singularity images are stored. 
When using a computing cluster it must be a shared folder accessible from all computing nodes.
```

* Execute the following for LSF on juno, setting the output directory to `test_output` (below) or another subdirectory

```
nextflow run pipeline.nf --somatic --germline --outDir test_output --mapping test_inputs/lsf/test_make_bam_and_qc.tsv --pairing test_inputs/lsf/test_make_bam_and_qc_pairing.tsv -profile juno
```