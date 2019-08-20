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
└── alignment_qc.txt
```

These outputs are:
- `fastp` (folder): An HTML report for each FASTQ file pair per sample.
- `alfred` (folder): A per-sample and per-readgroup BAM file alignment metrics in text and PDF files.
- `collectshsmetrics` (folder): For exomes, per-sample hybridisation-selection metrics in the.
- `conpair` (folder): Per tumor-normal-pair contamination and sample concordance estimates.
- `alignment_qc.txt`: Aggregated read-alignments statistics file, from the `alfred` and `collectshsmetrics` folders.

## Somatic data

The result of the somatic analyses is output in summarized forms in the `somatic` folder: 

```shell
somatic
├── facets
├── mutsig
├── merged_all_neoantigen_predictions.txt
├── merged_armlevel.tsv
├── merged_genelevel_TSG_ManualReview.txt
├── merged_hisens.cncf.txt
├── merged_hisensPurity_out.txt
├── merged_hisens.seg
├── merged.maf
├── merged_metadata.tsv
├── merged_purity.cncf.txt
├── merged_purity.seg
└── merged.vcf.gz
```

These outputs are:
- `facets` (folder): Individual copy-number profiles from FACETS, per tumor-normal pair.
- `mutsig` (folder): Individual mutational signature decomposition per tumor-normal pair.
- `merged_all_neoantigen_predictions.txt`: Neoantigen predictions from NetMHCpan for all samples.
- `merged_armlevel.tsv`, `merged_genelevel_TSG_ManualReview.txt`, `merged_hisens.cncf.txt`, `merged_hisensPurity_out.txt`, `merged_purity.cncf.txt`, and `merged_purity.seg`, summarized arm- and gene-level output from Facets, as well as FACETS- and IGV-style segmentation files.
- `merged.maf`: Filtered mutations from MuTect2 and Strelka2, annotated with mutational effects, neoantigen predictions, and zygosity, as [described elsewhere](variant-annotation-and-filtering.md#somatic-snvs-and-indels).
- `merged.vcf.gz`: All structural variants detected by Delly and Manta.

## Germline data

The result of the somatic analyses is output in summarized forms in the `germline` folder: 

```shell
germline/
├── merged.maf
└── merged.vcf.gz
```

These outputs are:
- `merged.maf`: Filtered mutations from HaplotypeCaller and Strelka2, annotated with mutational effects and zygosity, as [described elsewhere](variant-annotation-and-filtering.md#germline-snvs-and-indels).
- `merged.vcf.gz`: All structural variants Delly and Manta.

## Extended Outputs

TBD
