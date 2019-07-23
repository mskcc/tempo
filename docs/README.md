---
home: true
heroImage: /vaporware_logo.jpg
actionText: Get Started →
actionLink: /installation/
features:
- title: Simplicity First
  details: Minimal setup with markdown-centered project structure helps you focus on writing.
- title: Simplicity First
  details: Minimal setup with markdown-centered project structure helps you focus on writing.
- title: Simplicity First
  details: Minimal setup with markdown-centered project structure helps you focus on writing.
footer: MIT Licensed | Copyright © 2019-present
---

# Vaporware documentation

Vaporware is a computational pipeline for processing data from paired-end exome and genome sequencing of human cancer samples with matched normals. Its components are containerized and the pipeline runs on the [Juno high-performance computing cluster](http://hpc.mskcc.org/) at Memorial Sloan Kettering Cancer Cencer and on [Amazon Web Services (AWS)](https://aws.amazon.com). The pipeline was written by members of the [Center for Molecular Oncology](https://cmo.mskcc.org).

These pages contain instructions on how to run the Vaporware pipeline. It also contains documentation on the bioinformatic components in the pipeline, some motivation for various parameter choices, plus sources and processing of all reference materials use. 

If there are any questions or comments, you are welcome to [raise an issue](https://github.com/mskcc/vaporware/issues/new?title=[User%20question]).

<small>Note: Vaporware currently only supports human samples, has only been tested for exome and genome sequencing experiments, and all references files are in build GRCh37 of the human genome.</small>

---

## Table of contents

### 1. Getting started

#### 1.1. [Installation](installation.md)
* [Set up on Juno](juno-setup.md)
* [Set up on AWS](aws-setup.md)

#### 1.2. [Usage](usage.md)
* [Nextflow basics](nextflow-basics.md)
* [Working with containers](working-with-containers.md)
* [Run pipeline](run-pipeline.md)

### 2. Pipeline contents

#### 2.1. [Bioinformatic components](bioinformatic-components.md)
* [Read alignment](bioinformatic-components.md#read-alignment)
* [Somatic analyses](bioinformatic-components.md#somatic-analyses)
* [Germline analyses](bioinformatic-components.md#germline-analyses)

#### 2.2. [Reference files](reference-files.md)
* [Intervals](intervals.md)

#### 2.3. [Variant annotation and filtering](variant-annotation-and-filtering.md)
* [Panel of normals for exomes](wes-panel-of-normals.md)

### 3. Help and other resources
* [Troubleshooting](troubleshooting.md)
* [AWS glossary](aws-glossary.md)

### 4. Contributing
* [Contributing to Vaporware](contributing-to-vaporware.md)
---