# Outputs

All paths below are relative to the base directory `outDir` as described in the [run instructions](running-the-pipeline.md).
```shell
outDir
├── bams
├── qc
├── somatic
└── germline
```

## BAM Files 

The `bams` folder contains the final aligned and post-processed BAM files along with index files.

## QC Outputs

FASTQ file, read alignment and basic BAM file QC is in the `qc` directory:

```shell
qc
├── alfred
├── collecthsmetrics
├── conpair
├── fastp
├── alignment_qc.txt
├── concordance_qc.txt
└── contamination_qc.txt
```

These outputs are:
- `fastp` (folder): An HTML report for each FASTQ file pair per sample.
- `alfred` (folder): A per-sample and per-readgroup BAM file alignment metrics in text and PDF files.
- `collectshsmetrics` (folder): For exomes, per-sample hybridisation-selection metrics in the.
- `conpair` (folder): Per tumor-normal-pair, the Conpair-generated SNP pileup files.
- `alignment_qc.txt`: Aggregated read-alignments statistics file, from the `alfred` and `collectshsmetrics` folders.
- `concordance_qc.txt`: Aggregated tumor-normal sample concordance estimates from Conpair.
- `contamination_qc.txt`: Aggregated sample contamination estimates from Conpair.

## Somatic data

The result of the somatic analyses is output in summarized forms in the `somatic` folder: 

```shell
somatic
├── facets
├── mut_somatic.maf
├── mut_somatic_neoantigens.txt
├── cna_armlevel.txt
├── cna_genelevel.txt
├── cna_hisens.seg
├── cna_purity.seg
├── cna_facets_run_info.txt
├── sv_somatic.vcf.{gz,.gz.tbi}
└── sample_data.txt
```

These outputs are:
- `facets` (folder): Individual copy-number profiles from FACETS, per tumor-normal pair.
- `mut_somatic.maf`: Filtered mutations from MuTect2 and Strelka2, annotated with mutational effects, neoantigen predictions, and zygosity, as [described elsewhere](variant-annotation-and-filtering.md#somatic-snvs-and-indels).
- `mut_somatic_neoantigens.txt`: Neoantigen predictions from NetMHCpan for all samples.
- `cna_armlevel.txt`, `cna_genelevel.txt`, and `cna_hisens.seg`, `cna_purity.seg`, and `cna_facets_run_info.txt`, summarized arm- and gene-level output from Facets, as well as IGV-style segmentation files and Facets run information.
- `sv_somatic.vcf.gz`: All structural variants detected by Delly and Manta.
- `sample_data.txt` : Merged metadata across samples and analyses.

## Germline data

The result of the germline analyses is output in summarized forms in the `germline` folder: 

```shell
germline/
├── mut_germline.maf
└── sv_germline.{gz,.gz.tbi}
```

These outputs are:
- `mut_germline.maf`: Filtered mutations from HaplotypeCaller and Strelka2, annotated with mutational effects and zygosity, as [described elsewhere](variant-annotation-and-filtering.md#germline-snvs-and-indels).
- `sv_germline.vcf.gz`: All structural variants Delly and Manta.

## Extended Outputs

When run with the flag `--publishAll`, the pipeline will output additional intermediate data from select processes. These are:

```shell
somatic
├── mutations
    └── mutect2
    └── strelka2
├── structural_variants
    └── delly
    └── manta    
├── facets
└── lohhla

germline
├── mutations
    └── haplotypecaller
    └── strelka2
└── structural_variants
    └── delly
    └── manta
```

The `mutations` subdirectory contain VCFs with the unfiltered variant calls from the somatic and germline SNV/indel and SV callers. Additionally, these directories contain per-sample unfiltered MAF files generated in the `SomaticAnnotateMaf` and `GermlineAnnotateMaf` process, respectively. The `facets` subdirectory will contain the full arm- and gene-level outputs per sample. In the `lohhla` subdirectory the full LOHHLA LOH output metrics will be together with a PDF file with graphical output.

