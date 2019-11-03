---
home: true
heroImage: /tempoLogo.jpg
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

# Time-Efficient Mutational Profiling in Oncology (Tempo)

Tempo is a computational pipeline for processing data of paired-end whole-exome (WES) and whole-genome sequencing (WGS) of human cancer samples with matched normals. Its components are containerized and the pipeline runs on the [Juno high-performance computing cluster](http://mskcchpc.org/display/CLUS/Juno+Cluster+Guide) at Memorial Sloan Kettering Cancer Center and on [Amazon Web Services (AWS)](https://aws.amazon.com). The pipeline was written by members of the [Center for Molecular Oncology](https://www.mskcc.org/research-programs/molecular-oncology).

These pages contain instructions on how to run the Tempo pipeline. It also contains documentation on the bioinformatic components in the pipeline, some motivation for various parameter choices, plus an outline describing the reference resources used. 

If there are any questions or comments, you are welcome to [raise an issue](https://github.com/mskcc/tempo/issues/new?title=[User%20question]).

<small>Note: Tempo currently only supports human samples. The pipeline has only been tested for exome and genome sequencing experiments, and all reference files are in build GRCh37 of the human genome.</small>

---

## Table of Contents

### 1. Getting Started

#### 1.1. Setup
* [Installation](installation.md)
* [Setup on Juno](juno-setup.md)
* [Setup on AWS](aws-setup.md)

#### 1.2. Usage
* [Running the Pipeline](running-the-pipeline.md)
    * [Modifying or Resuming Pipeline Run](running-the-pipeline.md#modifying-or-resuming-pipeline-run)
    * [After Successful Run](running-the-pipeline.md#after-successful-run)
* [Nextflow Basics](nextflow-basics.md)
* [Working With Containers](working-with-containers.md)

#### 1.3 Outputs
* [BAM Files](outputs.md#bam-files)
* [QC Outputs](outputs.md#qc-outputs)
* [Somatic Data](outputs.md#somatic-data)
* [Germline Data](outputs.md#germline-data)
* [Extended Outputs](outputs.md#extended-outputs)

### 2. Pipeline contents

#### 2.1. Bioinformatic Components
* [Read Alignment](bioinformatic-components.md#read-alignment)
* [Somatic Analyses](bioinformatic-components.md#somatic-analyses)
* [Germline Analyses](bioinformatic-components.md#germline-analyses)
* [Quality Control](bioinformatic-components.md#quality-control)

#### 2.2. Reference Resources
* [Genome Assembly](reference-files.md#genome-assembly)
* [Genomic Intervals](reference-files.md#genomic-intervals)
* [RepeatMasker and Mappability Blacklist](reference-files.md#repeatmasker-and-mappability-blacklist)
* [Preferred Transcript Isoforms](reference-files.md#preferred-transcript-isoforms)
* [Hotspot Annotation](reference-files.md#hotspot-annotation.md)
* [OncoKB Annotation](reference-files.md#oncokb.md)
* [gnomAD](gnomad.md)
* [Panel of Normals for Exomes](wes-panel-of-normals.md)
* [Gene Annotation for SVs](gene-annotation-for-svs.md)

#### 2.3. Variant Annotation and Filtering
* [Somatic SNVs and Indels](variant-annotation-and-filtering.md#somatic-snvs-and-indels)
* [Germline SNVs and Indels](variant-annotation-and-filtering.md#germline-snvs-and-indels)
* [Somatic and Germline SVs](variant-annotation-and-filtering.md#somatic-and-germline-svs)

### 3. Help and Other Resources
* [Troubleshooting](troubleshooting.md)
* [AWS Glossary](aws-glossary.md)

### 4. Contributing
* [Contributing to Tempo](contributing-to-tempo.md)

### 5. Acknowledgements
* [Acknowledgements](acknowledgements.md)

<p align="center">
  <img src="./brandenburg5_allegro.jpg">
</p>

---
