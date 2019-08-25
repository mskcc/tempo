# Reference Resources

This and associated pages in this section provide details on provenance and generation of all reference files used in `pipeline.nf`. Usage of these files is defined in the [references configuration file](../conf/references.config).

::: tip Note
All reference files described herein are in assembly GRCh37/hg19 of the human genome.
:::

## Genome Assembly

Part of the [GATK bundle](https://software.broadinstitute.org/gatk/download/bundle), also available [here](https://console.cloud.google.com/storage/browser/gatk-legacy-bundles/b37). Vaporware uses the **human_g1k_v37_decoy** assembly of the genome.

## Genomic Intervals

BED files that specify the regions of the genome to consider for variant calling are specified in the [input files](running-the-pipeline.md#input-files).

### Exome Capture Platforms
For exomes, use BED file corresponding to the platform used for target capture. Currently, Vaporware supports:
- __AgilentExon_51MB__: SureSelectXT Human All Exon V4 from Agilent.
- __IDT_Exome__: xGen Exome Research Panel v1.0 from IDT.

::: tip Note
Contact us if you are interested in support for other sequencing assays or capture kits.
:::

The bait and target files are provided by the kit manufacturer. These are used to estimate bait- and target-level coverage metrics as well as for variant calling.

We add 5 bp to each end of exons in the target file to make sure splice site mutations can be called:
``` shell
bedtools slop \
    -g b37.chrom.sizes \
    -i targets.bed \
    -r 5 \
    -l 5 \
    > targets.plus5bp.bed
```

### Callable Regions for Genomes
For genomes, a list of "callable" regions from GATK's bundle is used. This is converted from an interval list to a BED file:
```shell
gatk IntervalListToBed \
    --INPUT b37_wgs_calling_regions.v1.interval_list \
    --OUTPUT b37_wgs_calling_regions.v1.bed
```

## RepeatMasker and Mappability Blacklist
BED files with genomic repeat and mappability information are used to annotate the VCFs with somatic and germline SNV/indels. These data are from [RepeatMasker](http://www.repeatmasker.org/) and the [ENCODE consortium](http://rohsdb.cmb.usc.edu/GBshape/ENCODE/index.html), and the files are retrieved from the [UCSC Genome Browser](https://genome.ucsc.edu) and parsed as such:

``` shell
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/rmsk.txt.gz
gunzip rmsk.txt.gz
cut -f6-8,12 rmsk.txt | \
    grep -e "Low_complexity" -e "Simple_repeat" | \
    sed 's/^chr//g'> rmsk_mod.bed
bgzip rmsk_mod.bed
tabix --preset bed rmsk_mod.bed.gz

wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/wgEncodeDacMapabilityConsensusExcludable.bed.gz
gunzip wgEncodeDacMapabilityConsensusExcludable.bed.gz
sed -i 's/^chr//g' wgEncodeDacMapabilityConsensusExcludable.bed
bgzip wgEncodeDacMapabilityConsensusExcludable.bed
tabix --preset bed wgEncodeDacMapabilityConsensusExcludable.bed.gz
``` 

## Preferred Transcript Isoforms
The `--custom-enst` argument to vcf2maf takes a list of preferred gene transcript isoforms which mutations are mapped onto. We supply a consensus list of [`isoform_overrides_at_mskcc` and `isoform_overrides_uniprot`](https://github.com/mskcc/vcf2maf/tree/master/data), generated as such:
``` r
t1 = readr::read_tsv('isoform_overrides_at_mskcc')
t2 = readr::read_tsv('isoform_overrides_uniprot')
t2 %>%
    dplyr::filter(gene_name %nin% t1$gene_name) %>%
    dplyr::bind_rows(., t1) %>%
    readr::write_tsv('isoforms')
```

## Hotspot Annotation
Three types of mutation hotspots are annotated in the somatic MAF. These include SNV, indel in linear space as well as SNV hotspots in 3D space. These are annotated with the [annotateMaf package](https://github.com/taylor-lab/annotateMaf). 

## OncoKB Annotation
Functional mutation effects and predicted oncogenicity of variants, as well as level of clinical actionability are from [OncoKB](https://oncokb.org) and annotated using the [OncoKB annotator](https://github.com/oncokb/oncokb-annotator).

## BRCA Exchange Annotation
Annotation of germline variants in _BRCA1_ and _BRCA2_ is carried out with the [annotateMaf package](https://github.com/taylor-lab/annotateMaf). This includes variant-level annotation from the ENIGMA consortium and ClinVar.

## Structural Variant Calling
Delly provides and takes as an argument a [file of regions](https://github.com/dellytools/delly/tree/master/excludeTemplates) to _exclude_ from variant calling. This excludes telomeres and centromeres from auto- and allosomes as well as any other contig.

For Manta, subtract these regions from a bed file of the whole genome to generate a list of regions to _include_. First clean up the file provided by Delly, since it is not in `bed` format:
``` shell
grep -Ev "chr|MT|GL00|NC|hs37d5" human.hg19.excl.tsv > human.hg19.excl.clean.bed
bedtools subtract -a b37.bed -b human.hg19.excl.clean.bed > b37.minusDellyExclude.bed
```
