# Bioinformatic Components

The 4 main sections of the pipeline are:

1. Read alignment
2. Somatic variant detection
3. Germline variant detection
4. Quality Ccontrol

Additionally, various QC metrics are generated. Below are described the separate modules tools used. The following diagram outlines the workflow:

<img id="diagram" src="./pipeline-flowchart.png"/>
<img id="dag" src="./dag.png"/>

<small>Note: The pipeline can be run with already-aligned BAM files as input, which avoids the first of these three modules.</small>

## Read Alignment

Tempo accepts as input sequencing reads from one or multiple FASTQ file pairs (corresponding to separate sequencing lanes) per sample, as [described](running-the-pipeline.md#the-mapping-file). These are aligned against the human genome reference using common practices, which include:

- **Alignment** using [BWA mem](http://bio-bwa.sourceforge.net/), followed by conversion to BAM file format and sorting using [samtools](https://samtools.github.io).
- **Merging** of BAM files across sequencing lanes using [samtools](https://samtools.github.io).
- **PCR-duplicate marking** using [GATK MarkDuplicates](https://software.broadinstitute.org/gatk).
- **Base-quality score recalibration** with [GATK BaseRecalibrator and ApplyBQSR](https://software.broadinstitute.org/gatk/).

## Somatic Analyses

- **SNVs and indels** are called using [MuTect2](https://software.broadinstitute.org/gatk/documentation/tooldocs/4.beta.4/org_broadinstitute_hellbender_tools_walkers_mutect_Mutect2.php) and [Strelka2](https://github.com/Illumina/strelka). Subsequently, they are combined, annotated and filtered as described [in the section on variant annotation and filtering](variant-annotation-and-filtering.md#somatic-snvs-and-indels).
- **Structural variants** are detected by multiple callers and then merged, filtered and annotated as described [in the section on variant annotation and filtering](variant-annotation-and-filtering.md#somatic-and-germline-svs). Somatic and germline variants are generated with [Delly](https://github.com/dellytools/delly), [Manta](https://github.com/Illumina/manta) and [SvABA](https://github.com/walaj/svaba). Somatic whole genome variants are generated with the same callers, along with [BRASS](https://github.com/cancerit/BRASS).
- **Copy-number analysis** is performed with [FACETS](https://github.com/mskcc/facets) and processed using [facets-suite](https://github.com/mskcc/facets-suite). Locus-specific copy-number, purity and ploidy estimates are integrated with the SNV/indel calls to perform clonality and zygosity analyses.
- **Microsatellite instability** is detected using [MSIsensor](https://github.com/ding-lab/msisensor).
- **HLA genotyping** is performed with [POLYSOLVER](https://software.broadinstitute.org/cancer/cga/polysolver).
- **LOH at HLA loci** is assessed with [LOHHLA](https://github.com/mskcc/lohhla).
- **Mutational signatures** are inferred with [https://github.com/mskcc/tempoSig](https://github.com/mskcc/tempoSig).
- **Neoantigen prediction** using estimates of class I MHC binding affinity is performed with [NetMHC 4.0](https://www.ncbi.nlm.nih.gov/pubmed/28978689) and integrated into the set of SNV/indel calls using [https://github.com/taylor-lab/neoantigen-dev](https://github.com/taylor-lab/neoantigen-dev) (_Note: this repository is currently private_).
- **Clonality using Structural Variants** is assessed using [SVclone](https://github.com/mcmero/SVclone). Joint calling is performed on both SNPs and structural variants, providing an alternative to the SNP clonality inference offered by Facets, while also assigning structural variants in the same clonal structure.
- **Structural variant signatures** are inferred with [signature.tools.lib](https://github.com/Nik-Zainal-Group/signature.tools.lib)
- **Homologous recombination deficiency** is assessed using [HRDetect from the signature.tools.lib package](https://github.com/Nik-Zainal-Group/signature.tools.lib)

## Germline Analyses

- **SNVs and indels** are called using [HaplotypeCaller](https://software.broadinstitute.org/gatk/documentation/tooldocs/4.0.8.0/org_broadinstitute_hellbender_tools_walkers_haplotypecaller_HaplotypeCaller.php) and [Strelka2](https://github.com/Illumina/strelka). Subsequently, they are combined, annotated and filtered as described [in the section on variant annotation and filtering](variant-annotation-and-filtering.md#germline-snvs-and-indels).
- **Structural variants** are detected by [Delly](https://github.com/dellytools/delly) and [Manta](https://github.com/Illumina/manta) then combined, filtered and annotated as described [in the section on variant annotation and filtering](variant-annotation-and-filtering.md#somatic-and-germline-svs).

## Quality Control

- **FASTQ QC metrics** are generated using [fastp](https://github.com/OpenGene/fastp).
- **BAM file QC metrics** are generated using [Alfred](https://github.com/tobiasrausch/alfred).
- **Hybridisation-selection metrics** are generated using [CollectHsMetrics](https://software.broadinstitute.org/gatk/documentation/tooldocs/4.beta.6/picard_analysis_directed_CollectHsMetrics.php). Only for exomes.
- **Contamination and concordance metrics** for tumor-normal pairs using [Conpair](https://github.com/mskcc/Conpair).

### MultiQC Report

A combined MultiQC report is produced at the BAM level, somatic level and cohort level of analysis, highlighting QC metrics and high-level summaries produced by Tempo. Documentation for the QC tool can be found [here](https://multiqc.info/docs/).

The Tempo team has internally derived `pass`/`warn`/`fail` thresholds to present in the QC report. This should be used with discretion by the analyst, as a failure in the QC report may not be a true failure, and so on. By default the following metrics are assessed and a `Status` value is produced in the report:

- **Tumor_Contamination** : Estimates percentage of cross-individual contamination in the tumor sample. Cross-contamination may be higher if the sample comes from the recipient of allograft tissue.
- **Normal_Contamination** : Estimates percentage of cross-individual contamination in the normal sample. Cross-contamination may be higher if the sample comes from the recipient of allograft tissue.
- **Concordance** : measures likelihood of samples coming from the same individual. Low values may indicate sample swap or contamination.
- **facets_qc** and **purity** : Report by [facets-preview](https://github.com/taylor-lab/facets-preview). A `facets_qc` value of `FALSE` indicates that the pair has failed Facets QC.
- **Fold Enrichment** : The fold by which the baited region has been amplified above genomic background. A low metric could indicate inefficiency of the bait selection kit during sample preparation.
- ** Target Bases 50X** : The fraction of all target bases achieving 50X or greater coverage. A low metric could indicate insufficient coverage.
