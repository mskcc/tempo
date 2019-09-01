# Bioinformatic Overview


The three main functions of the pipeline are:
- Sequencing read alignment
- Somatic variant detection
- Germline variant detection

Below are described the separate Nextflow processes and associated with each module and the tools used. The following diagram outlines the workflow:

<img id="diagram" src="./pipeline-flowchart.png"/>


<small>Note: The pipeline can be run with already-aligned BAM files as input, which avoids the first of these three modules.</small>

## Read alignment
Vaporware accepts as input sequencing reads from one or multiple FASTQ file pairs (corresponding to separate sequencing lanes) per sample, as [described](run-pipeline.md#the-mapping-file). These are aligned against the human genome using common practices.
* __AlignReads__: Alignment using [BWA mem](http://bio-bwa.sourceforge.net/), followed by conversion to BAM file format and sorting using [samtools](https://samtools.github.io). Additionally, this generates FASTQ QC metrics using [fastp](https://github.com/OpenGene/fastp).
* __MergeBams__: Merging of BAM files across sequencing lanes using [samtools](https://samtools.github.io).
* __MarkDuplicates__: PCR-duplicate marking using [GATK MarkDuplicates](https://software.broadinstitute.org/gatk).
* __CreateRecalibrationTable__ and __RecalibrateBam__: Base-quality score recalibration with [GATK BaseRecalibrator and ApplyBQSR](https://software.broadinstitute.org/gatk/).
* __Alfred__: BAM file QC using [Alfred](https://github.com/tobiasrausch/alfred).

## Somatic analyses
* SNVs and InDels are called using [MuTect2](https://software.broadinstitute.org/gatk/documentation/tooldocs/4.beta.4/org_broadinstitute_hellbender_tools_walkers_mutect_Mutect2.php) and [Strelka2](https://github.com/Illumina/strelka).
* Structural variants are detected by [Delly](https://github.com/dellytools/delly) and [Manta](https://github.com/Illumina/manta).
* Copy-number analysis is performed with [FACETS](https://github.com/mskcc/facets) and processed using [facets-suite](https://github.com/mskcc/facets-suite).
* Microsatellite instability is detected using [MSIsensor](https://github.com/ding-lab/msisensor).
* HLA genotyping with [POLYSOVLER](https://software.broadinstitute.org/cancer/cga/polysolver).
* Allele-specific HLA LOH with [LOHHLA](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5720478/).
* Mutational signatures with https://github.com/mskcc/mutation-signatures
* Class I MHC Binding Affinity Prediction with [NetMHC 4.0](https://www.ncbi.nlm.nih.gov/pubmed/28978689)


## Germline analyses
* SNVs and InDels are called using [MuTect2](https://software.broadinstitute.org/gatk/documentation/tooldocs/4.beta.4/org_broadinstitute_hellbender_tools_walkers_mutect_Mutect2.php) and [Strelka2](https://github.com/Illumina/strelka).
* Structural variants are detected by [Delly](https://github.com/dellytools/delly) and [Manta](https://github.com/Illumina/manta).
* We currently do not have a special module dedicated to germline CNVs besides that provided by the germline SV calls. 


## WGS vs. WES

Contact us if you are interest in support for other sequencing assays or capture kits.