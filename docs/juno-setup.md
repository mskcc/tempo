# Juno setup

The Juno compute cluster is accessible to most researchers at MSKCC. If you do not have an account on Juno or have other questions about their services, contact [HPC](http://hpc.mskcc.org/contact-us).

## Java
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

## Singularity
As described in the page about [containers](working-with-containers.md), execution of Vaporware on Juno requires Singularity. 

#### Installation
We recommend using Singularity version 3.1.1, as such:
```shell
module load singularity/3.1.1
```
or:
```shell
export PATH=/opt/local/singularity/3.1.1/bin:$PATH
```
The call `which singularity` should return `/opt/local/singularity/3.1.1/bin/singularity` if you have done this correctly.

#### Containers
By default, Singularity will look for and download image files to `~/.singularity`. In order to save time and space, you can use image files stored in a common cache directory by setting the environment variable `NXF_SINGULARITY_CACHEDIR` to the directory `/juno/work/taylorlab/cmopipeline/sandbox/singularity_images` (TO BE CHANGED!). You can put this in your bash profile:
```shell
export NXF_SINGULARITY_CACHEDIR=/juno/work/taylorlab/cmopipeline/sandbox/singularity_images
```
If you want to maintain your own cache of images, set this to your directory of choice.


