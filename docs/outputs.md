# Outputs

All paths below are relative to the base directory `outDir` as described in the [run instructions](running-the-pipeline.md).
```shell
outDir
├── bams
├── cohort_level
├── germline
└── somatic
```

## BAM Files 

The `bams` folder contains the final aligned and post-processed BAM files along with index files.
It also contains FASTQ QC and basic BAM file QC.
```shell
outDir/bams/
├── DU874145-N
│   ├── DU874145-N.bam
│   ├── DU874145-N.bam.bai
│   ├── alfred
│   ├── collecthsmetrics
│   ├── fastp
│   └── pileup
├── DU874145-T
│   ├── DU874145-T.bam
│   └── DU874145-T.bam.bai
...
```

These outputs are:
- `fastp`: A HTML report for each FASTQ lane pair per sample.
- `alfred`: A per-sample and per-readgroup BAM file alignment metrics in text and PDF files.
- `collectshsmetrics`: For exomes, per-sample hybridisation-selection metrics in the.
- `pileup`: Per tumor-normal-pair, the Conpair-generated SNP pileup files.

## Somatic data

The result of the somatic analyses is output in summarized forms in the `somatic` folder: 

```shell
outDir/somatic
├── DU874145-T__DU874145-N
│   ├── combined_mutations
│   ├── combined_svs
│   ├── conpair
│   ├── delly
│   ├── facets
│   ├── lohhla
│   ├── manta
│   ├── meta_data
│   ├── mutect2
│   ├── neoantigen
│   └── strelka2
└── DU874146-T__DU874146-N
    ├── combined_mutations
    ├── combined_svs
    ├── delly
    ├── facets
    ├── lohhla
    ├── manta
    ├── meta_data
    ├── mutect2
    ├── neoantigen
    └── strelka2
```

These outputs are:
- `combined_mutatations`: unfiltered and final filtered maf per tumor-normal pair.
  - `*.somatic.unfiltered.maf`: Unfiltered mutations `generated in the SomaticAnnotateMaf`.
  - `*.somatic.final.maf`: Filtered mutations from MuTect2 and Strelka2, annotated with mutational effects, neoantigen predictions, and zygosity, as [described elsewhere](variant-annotation-and-filtering.md#somatic-snvs-and-indels).
- `combined_svs`: Combined Delly and Manta SV calls.
- `conpair`: Per tumor-normal-pair, the Conpair-generated concordance and contamination files.
- `delly`: Delly output.
- `facets`: Individual copy-number profiles from FACETS, per tumor-normal pair.
- `lohhla`: LOHHLA output.
- `manta`: Manta output.
- `meta_data`: Summarized meta_data file which includes the following results:
  - Purity and Ploidy
  - WGS Status
  - MSI information including MSI_Total_Sites, MSI_Somatic_Sites, MSIscore
  - Number of Mutations
  - All 30 Mutational Signatures
  - HLA genotyping
  - TMB
- `mutect2`: Manta output.
- `neoantigens`: Neoantigen predictions from NetMHCpan per sample.
- `strelka2`: Manta output.

::: warning Be aware
* LOHHLA is temporarily disabled due to a bug need future investigation. It will be enabled again in the future release.
:::

## Germline data

The result of the germline analyses is output in the `germline` folder:

```shell
outDir/germline/
├── DU874145-N
│   ├── combined_mutations
│   ├── combined_svs
│   ├── delly
│   ├── haplotypecaller
│   ├── manta
│   └── strelka2
└── DU874146-N
    ├── combined_mutations
    ├── combined_svs
    ├── delly
    ├── haplotypecaller
    ├── manta
    └── strelka2
```

These outputs are:
- `combined_mutatations`: unfiltered and final filtered maf per tumor-normal pair.
  - `*.germline.unfiltered.maf`: Unfiltered mutations `generated in the GermlineAnnotateMaf`.
  - `*.germline.final.maf`: Filtered mutations from HaplotypeCaller and Strelka2, annotated with mutational effects and zygosity, as [described elsewhere](variant-annotation-and-filtering.md#germline-snvs-and-indels).
- `combined_svs`: Combined Delly and Manta SV calls.
- `delly`: Delly output.
- `manta`: Manta output.
- `strelka2`: Manta output.

## Cohort Level Outputs

When run with the flag `--aggregate`, the pipeline will output aggregate all samples together for each processes and output as a single file for each processes. The files are:

```shell
outDir/cohort_level/
├── default_cohort
│   ├── alignment_qc.txt
│   ├── cna_armlevel.txt
│   ├── cna_facets_run_info.txt
│   ├── cna_genelevel.txt
│   ├── cna_hisens_run_segmentation.seg
│   ├── cna_purity_run_segmentation.seg
│   ├── concordance_qc.txt
│   ├── contamination_qc.txt
│   ├── DNA.IntegerCPN_CI.txt
│   ├── HLAlossPrediction_CI.txt
│   ├── mut_germline.maf
│   ├── mut_somatic.maf
│   ├── mut_somatic_neoantigens.txt
│   ├── sample_data.txt
│   ├── sv_germline.vcf.gz
│   ├── sv_germline.vcf.gz.tbi
│   ├── sv_somatic.vcf.gz
│   └── sv_somatic.vcf.gz.tbi
├── cohort2
│   ├── alignment_qc.txt
│   ├── cna_armlevel.txt
│   ├── cna_facets_run_info.txt
│   ├── cna_genelevel.txt
...
```

These outputs are just naively concatenated together from per sample output files (duplicated header are removed).
