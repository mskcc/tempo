# mskcc/tempo: Usage

## Introduction

The TEMPO (Tumor Exome aNalysis Pipeline Or-iented) pipeline is a comprehensive nf-core workflow for processing paired-end whole exome or genome sequencing data from tumor and normal samples. It performs quality control, alignment, somatic and germline variant calling, copy number analysis, microsatellite instability assessment, and HLA typing.

## Samplesheet Input

### Format and Requirements

The pipeline requires a CSV (comma-separated values) samplesheet specified via the `--input` parameter. This file defines all sample information and sequencing data paths required for the analysis.

### Columns

| Column | Description | Valid Values |
|--------|-------------|--------------|
| `patient` | Unique patient identifier | Alphanumeric string |
| `sample` | Unique sample identifier within patient | Alphanumeric string |
| `status` | Sample type classification | `0` (normal), `1` (tumor) |
| `sex` | Biological sex for sex chromosome analysis | `XX`, `XY`, `NA` |
| `lane` | Sequencing lane identifier | Alphanumeric string (e.g., `L001`, `L002`) |
| `fastq_1` | Path to first read FASTQ file | Full absolute path to `*_R1.fastq.gz` |
| `fastq_2` | Path to second read FASTQ file | Full absolute path to `*_R2.fastq.gz` |

### Example Samplesheet

```csv
patient,sample,status,sex,lane,fastq_1,fastq_2
patient1,sample1_normal,0,XX,L001,/path/to/normal_R1.fastq.gz,/path/to/normal_R2.fastq.gz
patient1,sample1_tumor,1,XX,L001,/path/to/tumor_R1.fastq.gz,/path/to/tumor_R2.fastq.gz
patient2,sample2_normal,0,XY,L001,/path/to/normal_R1.fastq.gz,/path/to/normal_R2.fastq.gz
patient2,sample2_tumor,1,XY,L001,/path/to/tumor_R1.fastq.gz,/path/to/tumor_R2.fastq.gz
patient2,sample2_tumor,1,XY,L002,/path/to/tumor_R1.fastq.gz,/path/to/tumor_R2.fastq.gz
```

### Validation Rules

- **One normal per patient**: Each patient MUST have exactly one normal sample (status=0)
- **One or more tumors per patient**: Each patient MUST have at least one tumor sample (status=1)
- **Multi-lane samples**: If a sample was sequenced across multiple lanes, provide separate rows with identical patient, sample, and status but different lane values and FASTQ paths. The pipeline will automatically concatenate reads from multiple lanes.
- **Path requirements**: FASTQ file paths must be absolute paths to gzip-compressed files (.fastq.gz or .fq.gz)
- **File existence**: All specified FASTQ files must exist and be readable before running the pipeline

### Sample Relationships

The pipeline processes tumor-normal pairs. Each tumor sample is paired with the single normal sample from the same patient for somatic variant calling and quality control. For patients with multiple tumor samples, each tumor is independently compared to the same normal sample.

## Running the Pipeline

### Basic Execution

The typical command for running the pipeline is as follows:

```bash
nextflow run mskcc/tempo --input samplesheet.csv --outdir results -profile singularity
```

This will launch the pipeline with the `singularity` configuration profile, which is recommended for HPC environments.

### On Juno Cluster

For users on the MSK Juno cluster with preconfigured settings:

```bash
nextflow run mskcc/tempo --input samplesheet.csv --outdir results -profile juno
```

### Test Run

For testing the pipeline with provided sample data:

```bash
nextflow run mskcc/tempo -profile test,docker --outdir results
```

The test profile includes pre-configured inputs and reference files for validation.

### Resume Previous Run

To resume an interrupted run from the last successful task:

```bash
nextflow run mskcc/tempo --input samplesheet.csv --outdir results -profile singularity -resume
```

### Working Directory Structure

The pipeline will create the following files in your working directory:

```bash
work                # Directory containing the nextflow working files
results             # Finished results in specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

### Common Command-Line Parameters

| Parameter | Description | Default | Required |
|-----------|-------------|---------|----------|
| `--input` | Path to samplesheet CSV | - | Yes |
| `--outdir` | Output directory | `./results` | No |
| `-profile` | Configuration profile (singularity, docker, juno, test) | - | Yes |
| `-resume` | Resume from last successful task | false | No |
| `-r` | Specific pipeline version | Latest | No |

### Updating the Pipeline

When you run the pipeline command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available, even if the pipeline has been updated. To ensure you're running the latest version:

```bash
nextflow pull mskcc/tempo
```

### Reproducibility

It is a best practice to specify a pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used, allowing for reproducible analyses.

To specify a version, first visit the [mskcc/tempo releases page](https://github.com/mskcc/tempo/releases) and find the desired version number (e.g., `1.3.1`). Then use the `-r` flag when running:

```bash
nextflow run mskcc/tempo -r 1.3.1 --input samplesheet.csv --outdir results -profile singularity
```

The version number will be logged in execution reports and MultiQC output for future reference.

## Reference Genomes

### Default Genome

GRCh37 (hg19) is the default reference genome used by the pipeline. All results are reported against this build unless otherwise specified.

### Custom Reference Genomes

To use a different reference genome, you must provide the following required files:

### Required Reference Files

| Parameter | Description | Format |
|-----------|-------------|--------|
| `--fasta` | Reference genome FASTA file | FASTA (.fa or .fasta) |
| `--fasta_fai` | FASTA index file | Generated with `samtools faidx` |
| `--dict` | Dictionary file | Generated with `picard CreateSequenceDictionary` |
| `--bwa_index` | BWA index files | Prefix for `*.amb`, `*.ann`, `*.bwt`, `*.pac`, `*.sa` files |
| `--dbsnp` | dbSNP known variants | VCF (.vcf.gz) |
| `--known_indels` | Known indel locations | VCF (.vcf.gz) (e.g., Mills and 1000G gold standard) |
| `--germline_resource` | Germline variants for contamination filtering | VCF (.vcf.gz) (e.g., gnomAD) |
| `--intervals` | Target regions for analysis | BED or interval_list format |

### TEMPO-Specific Reference Files

| Parameter | Description | Format |
|-----------|-------------|--------|
| `--facets_vcf` | Common SNP VCF for FACETS copy number analysis | VCF (.vcf.gz) |
| `--msi_sensor_list` | Microsatellite list for MSIsensor-pro scoring | List format |
| `--vep_cache` | VEP (Variant Effect Predictor) annotation cache | Directory |
| `--hla_fasta` | HLA reference sequences for HLA typing | FASTA |

### Example Reference Configuration

```bash
nextflow run mskcc/tempo \
  --input samplesheet.csv \
  --outdir results \
  -profile singularity \
  --fasta /path/to/GRCh37.fa \
  --fasta_fai /path/to/GRCh37.fa.fai \
  --dict /path/to/GRCh37.dict \
  --bwa_index /path/to/bwa_index/GRCh37 \
  --dbsnp /path/to/dbsnp_146.vcf.gz \
  --known_indels /path/to/Mills_and_1000G_gold_standard.indels.vcf.gz \
  --germline_resource /path/to/af-only-gnomad.vcf.gz \
  --intervals /path/to/exome.bed \
  --facets_vcf /path/to/facets_snps.vcf.gz \
  --msi_sensor_list /path/to/msi_sensor.list \
  --vep_cache /path/to/vep_cache \
  --hla_fasta /path/to/hla.fasta
```

## Pipeline Options

### Assay Type

Specify the sequencing assay type. This affects the analysis approach and reference regions used:

```bash
--assay_type exome    # (default) Whole exome sequencing
--assay_type genome   # Whole genome sequencing
```

### Skip Options

Skip specific analysis modules to reduce runtime or for debugging purposes:

```bash
--skip_somatic_snv      # Skip Mutect2 and Strelka2 somatic SNV calling
--skip_somatic_sv       # Skip Manta and Delly structural variant calling
--skip_facets           # Skip FACETS copy number analysis
--skip_msi              # Skip MSIsensor-pro microsatellite instability analysis
--skip_polysolver       # Skip Polysolver HLA typing
--skip_lohhla           # Skip LOHHLA HLA loss of heterozygosity analysis
--skip_germline_snv     # Skip HaplotypeCaller germline SNV calling
--skip_qc               # Skip Conpair concordance/contamination check
--skip_multiqc          # Skip MultiQC report generation
```

### Quality Control Options

```bash
--min_reads_unmapped        # Minimum percentage of unmapped reads to flag QC warning
--contamination_threshold   # Conpair contamination threshold for warning (default: 0.05)
```

### Variant Calling Options

```bash
--mutect2_extra_args        # Extra arguments to pass to Mutect2
--strelka_extra_args        # Extra arguments to pass to Strelka2
--manta_extra_args          # Extra arguments to pass to Manta
```

## Core Nextflow Arguments

> **NB:** These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen).

### `-profile`

Use this parameter to choose a configuration profile. Profiles provide configuration presets for different compute environments and container technologies.

Available profiles include:

- `test` - Complete configuration for automated testing with test data
- `docker` - Use Docker containers (recommended for local machines)
- `singularity` - Use Singularity containers (recommended for HPC clusters)
- `juno` - MSK Juno cluster specific configuration
- `podman` - Use Podman containers
- `shifter` - Use Shifter containers (NERSC)
- `charliecloud` - Use Charliecloud containers
- `conda` - Use Conda environment (not recommended for reproducibility)

Multiple profiles can be combined: `-profile test,docker` (order matters - later profiles override earlier ones).

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from pipeline steps where inputs are identical, continuing from the last successful task. This is useful for recovering from temporary failures.

```bash
nextflow run mskcc/tempo --input samplesheet.csv -resume
```

You can also resume a specific named run:

```bash
nextflow run mskcc/tempo --input samplesheet.csv -resume [run-name]
```

Use `nextflow log` to view previous run names.

### `-c`

Specify a custom Nextflow configuration file to override default settings:

```bash
nextflow run mskcc/tempo --input samplesheet.csv -c custom.config
```

See the [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information.

## Configuration and Advanced Options

### Nextflow Configuration File

Store pipeline parameters in a configuration file for consistent and reproducible runs across multiple executions:

```groovy
// nextflow.config
params {
  input = 'samplesheet.csv'
  outdir = 'results'
  assay_type = 'exome'

  // Reference files
  fasta = '/path/to/GRCh37.fa'
  fasta_fai = '/path/to/GRCh37.fa.fai'
  dict = '/path/to/GRCh37.dict'
  bwa_index = '/path/to/bwa_index/GRCh37'
  dbsnp = '/path/to/dbsnp_146.vcf.gz'
  known_indels = '/path/to/Mills_and_1000G_gold_standard.indels.vcf.gz'
  germline_resource = '/path/to/af-only-gnomad.vcf.gz'
  intervals = '/path/to/exome.bed'
  facets_vcf = '/path/to/facets_snps.vcf.gz'
  msi_sensor_list = '/path/to/msi_sensor.list'
  vep_cache = '/path/to/vep_cache'
  hla_fasta = '/path/to/hla.fasta'
}

process {
  executor = 'slurm'
  queue = 'default'
  memory = '4 GB'
  cpus = 4
}
```

Then run with:

```bash
nextflow run mskcc/tempo -profile singularity
```

### Resource Requests

Each step in the pipeline has default CPU, memory, and time requirements. If a process fails with specific error codes, Nextflow will automatically retry with increased resources (2x then 3x the original).

To globally adjust resources:

```bash
nextflow run mskcc/tempo --input samplesheet.csv --max_memory 200GB --max_cpus 16
```

### Updating Containers

To use a different version of a specific tool, create a custom configuration file:

```nextflow
process {
    withName: PROCESS_NAME {
        container = 'quay.io/biocontainers/tool:version'
    }
}
```

Then pass it to the pipeline:

```bash
nextflow run mskcc/tempo -c custom.config --input samplesheet.csv
```

### Running in the Background

Use Nextflow's background mode or terminal multiplexers to run long pipelines:

```bash
# Using Nextflow background mode
nextflow run mskcc/tempo --input samplesheet.csv -bg

# Using screen
screen -S tempo_run
nextflow run mskcc/tempo --input samplesheet.csv
# Detach with Ctrl+A then D
```

### Nextflow Memory Configuration

To limit Nextflow JVM memory usage, add to your shell profile (`~/.bashrc` or `~/.bash_profile`):

```bash
export NXF_OPTS='-Xms1g -Xmx4g'
```

## Troubleshooting

### Common Issues

- **Missing reference files**: Ensure all `--*` reference parameters are specified and point to valid, readable files
- **FASTQ path errors**: Verify that all paths in samplesheet are absolute paths to existing gzip-compressed files
- **Samplesheet validation failures**: Check that each patient has exactly one normal (status=0) and at least one tumor (status=1) sample
- **Out of memory errors (exit code 137)**: Increase memory allocation with `--max_memory`
- **Container issues**: Ensure Docker or Singularity is properly installed and configured
- **Resume failures**: Clean the `work/` directory if encountering resume issues, then restart from the beginning

### Getting Help

For additional support:
- Visit the [mskcc/tempo GitHub repository](https://github.com/mskcc/tempo)
- Check existing issues and discussions
- Review the [nf-core documentation](https://nf-co.re/)
- Consult the pipeline's troubleshooting guide
