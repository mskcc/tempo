# Bioinformatic components
The three main functions of the pipeline are:
- Sequencing read alignment
- Somatic variant detection
- Germline variant detection

<small>Note that the pipeline can be run with already aligned BAM files as input.</small>

## Read alignment
Vaporware accepts as input sequencing reads from one or multiple FASTQ file pairs (corresponding to separate sequencing lanes) per sample.
- Reads are aligned to the human genome using [BWA mem](http://bio-bwa.sourceforge.net/).
- Alignments are sorted, converted to BAM file format, and merged across sequencing lanes using [samtools](https://samtools.github.io).
- PCR duplicate marking and base-quality score recalibration are done with tools from the [GATK suite](https://software.broadinstitute.org/gatk/).

---
## Somatic analyses
- SNVs and indels are called using [MuTect2](https://software.broadinstitute.org/gatk/documentation/tooldocs/4.beta.4/org_broadinstitute_hellbender_tools_walkers_mutect_Mutect2.php) and [Strelka2](https://github.com/Illumina/strelka).
- Structural variants are detected by [Delly](https://github.com/dellytools/delly) and [Manta](https://github.com/Illumina/manta).
- Copy-number analysis is performed with [FACETS](https://github.com/mskcc/facets) and processed using [facets-suite](https://github.com/mskcc/facets-suite).
- Microsatellite instability is detected using [MSIsensor](https://github.com/ding-lab/msisensor).

---
## Germline analyses