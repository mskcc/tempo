# mskcc/tempo

[![Nextflow](https://img.shields.io/badge/nextflow%20-%E2%89%A524.04.2-23aa62.svg?colorB=0c3173)](https://www.nextflow.io/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/singularity/)
[![run with apptainer](https://img.shields.io/badge/run%20with-apptainer-1d355c.svg?labelColor=000000)](https://apptainer.org/)

## Introduction

**Tempo** (Time-Efficient Mutational Profiling in Oncology) is a comprehensive Nextflow pipeline for processing whole-exome and whole-genome sequencing (WES/WGS) data from tumor-normal pairs in cancer genomics research. Developed by the MSKCC Center for Molecular Oncology, Tempo implements nf-core best practices and supports both high-performance computing environments (Juno HPC at MSKCC with SLURM + Singularity) and cloud deployment (AWS).

The pipeline performs end-to-end analysis including quality control, read alignment, somatic variant calling, structural variant detection, copy number analysis, MSI detection, HLA typing, and comprehensive reporting.

## Pipeline Summary

Tempo performs the following analysis steps:

1. **Read Quality Control** — FastQC quality assessment of raw sequencing reads
2. **Adapter Trimming** — Adapter and low-quality base removal using fastp
3. **Read Alignment** — Alignment to reference genome using BWA-MEM2
4. **Post-Alignment Processing** — Duplicate marking (MarkDuplicates) and base quality score recalibration (BQSR)
5. **Somatic SNV/Indel Calling** — Variant calling using Mutect2 and Strelka2
6. **Somatic Structural Variant Calling** — SV detection using Manta and Delly
7. **Copy Number Analysis** — Copy number segment analysis using FACETS
8. **Microsatellite Instability (MSI) Detection** — MSI status determination using MSIsensor-pro
9. **HLA Typing** — HLA allele inference using Polysolver
10. **HLA Loss of Heterozygosity (LOH)** — HLA-specific LOH analysis using LOHHLA
11. **Germline Variant Calling** — Germline variant discovery using HaplotypeCaller
12. **Quality Control and Reporting** — Sample concordance assessment (Conpair) and comprehensive QC report generation (MultiQC)

## Quick Start

1. **Install Nextflow** (version 24.04.2 or later)

   ```bash
   curl -s https://get.nextflow.io | bash
   ```

2. **Pull the pipeline**

   ```bash
   nextflow pull mskcc/tempo
   ```

3. **Run with test profile**

   ```bash
   nextflow run mskcc/tempo -profile test,docker
   ```

   Replace `docker` with `singularity` or `apptainer` as needed for your environment.

## Usage

For detailed usage instructions, parameters, and configuration options, see [docs/usage.md](docs/usage.md).

## Pipeline Output

Comprehensive documentation of output files and directory structure is available in [docs/output.md](docs/output.md).

## Reference Genomes

Tempo supports the following reference genomes:

- **GRCh37** (default, primary)
- **GRCh38** (supported)

## Pipeline DAG

The pipeline workflow can be visualized using the metro map visualization:

![Tempo Pipeline DAG](docs/images/tempo_metro_map.png)

## Credits

Tempo was developed and is maintained by the **MSKCC Center for Molecular Oncology**.

## Citations

For citations and references, see [CITATIONS.md](CITATIONS.md).

## License

This project is licensed under the MIT License.
