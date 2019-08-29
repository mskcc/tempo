---
home: true
heroImage: /vaporwareLogo.jpg
actionText: Get Started →
actionLink: /installation/
features:
- title: Reproducible Workflows
  details: Containerized workflows with Docker and Singularity
- title: Portable 
  details: Tailored for LSF and AWS 
- title: User-Friendly 
  details: Written to be quickly run and used by anyone in the CMO
footer: MIT Licensed | Copyright © 2019-present
---

# Computational Pipeline for Whole-Genome and Whole-Exome Sequencing

Tempo is a computational pipeline for processing data of paired-end whole-exome (WES) and whole-genome sequencing (WGS) of human cancer samples with matched normals. Its components are containerized and the pipeline runs on the [Juno high-performance computing cluster](http://hpc.mskcc.org/) at Memorial Sloan Kettering Cancer Center and on [Amazon Web Services (AWS)](https://aws.amazon.com). The pipeline was written by members of the [Center for Molecular Oncology](https://www.mskcc.org/research-programs/molecular-oncology).

These pages contain instructions on how to run the Tempo pipeline. It also contains documentation on the bioinformatic components in the pipeline, some motivation for various parameter choices, plus sources and processing of reference resources used. 

If there are any questions or comments, you are welcome to [raise an issue](https://github.com/mskcc/tempo/issues/new?title=[User%20question]).

<small>Note: Tempo currently only supports human samples, has only been tested for exome and genome sequencing experiments, and all references files are in build GRCh37 of the human genome.</small>

---

## Table of Contents

### 1. Getting Started

#### 1.1. [Installation](installation.md)
* [Set up on Juno](juno-setup.md)
* [Set up on AWS](aws-setup.md)

#### 1.2. Usage
* [Running the Pipeline](running-the-pipeline.md)
    * [Modifying or Resuming Pipeline Run](running-the-pipeline.md#modifying-or-resuming-pipeline-run)
    * [After Successful Run](running-the-pipeline.md#after-successful-run)
* [Nextflow Basics](nextflow-basics.md)
* [Working With Containers](working-with-containers.md)

#### 1.3 [Outputs](outputs.md)
* [BAM Files](outputs.md#bam-files)
* [QC Outputs](outputs.md#qc-outputs)
* [Somatic Data](outputs.md#somatic-data)
* [Germline Data](outputs.md#germline-data)
* [Extended Outputs](outputs.md#extended-outputs)

### 2. Pipeline contents

#### 2.1. [Bioinformatics Components](bioinformatics-components.md)
* [Read Alignment](bioinformatic-components.md#read-alignment)
* [Somatic Analyses](bioinformatic-components.md#somatic-analyses)
* [Germline Analyses](bioinformatic-components.md#germline-analyses)
* [Quality Control](bioinformatic-components.md#quality-control)

#### 2.2. [Reference Resources](reference-resources.md)
* [Genome Assembly](reference-resources.md#genome-assembly)
* [Genomic Intervals](reference-resources.md#genomic-intervals)
* [RepeatMasker and Mappability Blacklist](reference-resources.md#repeatmasker-and-mappability-blacklist)
* [Preferred Transcript Isoforms](reference-resources.md#preferred-transcript-isoforms)
* [Hotspot Annotation](reference-resources.md#hotspot-annotation.md)
* [OncoKB Annotation](reference-resources.md#oncokb.md)
* [gnomAD](gnomad.md)
* [Panel of Normals for Exomes](wes-panel-of-normals.md)

#### 2.3. [Variant Annotation and Filtering](variant-annotation-and-filtering.md)
* [Somatic SNVs and Indels](variant-annotation-and-filtering.md#somatic-snvs-and-indels)
* [Germline SNVs and Indels](variant-annotation-and-filtering.md#germline-snvs-and-indels)
* [Somatic and Germline SVs](variant-annotation-and-filtering.md#somatic-and-germline-svs)

### 3. Help and Other Resources
* [Troubleshooting](troubleshooting.md)
* [AWS Glossary](aws-glossary.md)

### 4. Contributing
* [Contributing to Tempo](contributing-to-tempo.md)
---