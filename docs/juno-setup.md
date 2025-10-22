# Juno Setup

The [Juno compute cluster](http://mskcchpc.org/display/CLUS/Juno+Cluster+Guide) is accessible to researchers within the CMO. If you do not have an account on Juno or have other questions about their services, contact [HPC](http://hpc.mskcc.org/contact-us). Juno uses the LSF job scheduler which Tempo is configured to work with.

## Temporary Files

The pipeline processes will write temporary files to the directory defined by the `TMPDIR` variable in the user environment on the compute node where it is active. On Juno, this should point to `/scratch`, so make sure you set this in your bash profile as such:
```shell
export TMPDIR=/scratch/username
```

::: danger Warning
Each compute node has a `/scratch` directory. If you inadvertently fill up this directory on a given node, processes that require a lot of space might fail on this node. If you supsect this to be the case, you can check your disk usage by doing `ssh -A nodename "du -hs /scratch/username"`. If necessary, you can clean up this directory by doing `ssh -A nodename "rm -rf /scratch/username/*"`.
:::

## Singularity Containers

As described in the page about [containers](working-with-containers.md), execution of Tempo on Juno requires Singularity. 

In order to save time and space, you can use image files stored in a common cache directory by setting the environment variable `NXF_SINGULARITY_CACHEDIR` to the directory `/juno/work/taylorlab/cmopipeline/singularity_images`. You can put this in your bash profile:

```shell
export NXF_SINGULARITY_CACHEDIR=/juno/work/taylorlab/cmopipeline/singularity_images
```

If you want to maintain your own cache of images, set this to your directory of choice, and pull/build the images. 

We recommend using Singularity version 3.1.1, as such:
```shell
module load singularity/3.1.1
```
or:
```shell
export PATH=/opt/local/singularity/3.1.1/bin:$PATH
```
The command `which singularity` should return `/opt/local/singularity/3.1.1/bin/singularity` if you have done this correctly. 

## Java Version

Nextflow requires Java version 8 or later. On Juno, you can load it using `module`:
```shell
module load java/jdk1.8.0_202
```
or put it in your `PATH` by inserting this into your bash profile:
```shell
export JAVA_HOME=/opt/common/CentOS_7/java/jdk1.8.0_202/
export PATH=$JAVA_HOME/bin:$PATH
```
The call `which java` should return `/opt/common/CentOS_7/java/jdk1.8.0_202/bin/java` if you have done this correctly.

## Test Your Environment

You can run a the pipeline on small test files to ensure that you are ready to run real data. If you experience any issues, something in your environment might be the reason. The following should take approximately 30 minutes and all tasks should succeed at first attempt:

```shell
nextflow run dsl2.nf \
    --mapping test_inputs/local/full_test_mapping.tsv \
    --pairing test_inputs/local/full_test_pairing.tsv \
    -profile test_singularity \
    --outDir results
    --aggregate
```
