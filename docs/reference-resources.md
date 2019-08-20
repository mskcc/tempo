# Reference resources

## RepeatMasker and mappability blacklist
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

## Preferred transcript isoforms
The `--custom-enst` argument to vcf2maf takes a list of preferred gene transcript isoforms which mutations are mapped onto. We supply a consensus list of [`isoform_overrides_at_mskcc` and `isoform_overrides_uniprot`](https://github.com/mskcc/vcf2maf/tree/master/data), generated as such:
``` r
t1 = readr::read_tsv('isoform_overrides_at_mskcc')
t2 = readr::read_tsv('isoform_overrides_uniprot')
t2 %>%
    dplyr::filter(gene_name %nin% t1$gene_name) %>%
    dplyr::bind_rows(., t1) %>%
    readr::write_tsv('isoforms')
```

## Hotspot annotation
Three types of mutation hotspots are annotated in the somatic MAF. These include SNV, indel in linear space as well as SNV hotspots in 3D space. These are annotated with the [annotateMaf package](https://github.com/taylor-lab/annotateMaf). 

## OncoKB
Functional mutation effects and predicted oncogenicity of variants, as well as level of clinical actionability are from [OncoKB](https://oncokb.org) and annotated using the [OncoKB annotator](https://github.com/oncokb/oncokb-annotator).


## GRCh37

### Genome
[GATK bundle](https://software.broadinstitute.org/gatk/download/bundle), also available [here](https://console.cloud.google.com/storage/browser/gatk-legacy-bundles/b37).

### SNV and indel calling
For genomes, use [`bed` file of "callable"](https://github.mskcc/vaporware/blob/master/docs/INTERVALS.md#genome) regions from GATK's bundle.

For exomes, use `bed` file corresponding to the platform used for target capture, see [documentation on intervals](https://github.mskcc/vaporware/blob/master/docs/INTERVALS.md#exome-capture-platform).


### Structural Variant (SV) calling
Delly provides and takes as an argument a [file of regions](https://github.com/dellytools/delly/tree/master/excludeTemplates) to _exclude_ from variant calling. This excludes telomeres and centromeres from auto- and allosomes as well as any other contig.

For Manta, subtract these regions from a bed file of the whole genome to generate a list of regions to _include_. First clean up the file provided by Delly, since it is not in `bed` format:
``` shell
grep -Ev "chr|MT|GL00|NC|hs37d5" human.hg19.excl.tsv > human.hg19.excl.clean.bed
bedtools subtract -a b37.bed -b human.hg19.excl.clean.bed > b37.minusDellyExclude.bed
```

## GRCh38
[GATK bundle](https://software.broadinstitute.org/gatk/download/bundle), also available [here](https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38).

**Note:** Support for hg38 is currently somewhat limited. Please raise an issue at https://github.com/mskcc/vaporware/issues if you would like to process data with GRch38. However, please note that hg38 reference files are easily available from UCSC. 

