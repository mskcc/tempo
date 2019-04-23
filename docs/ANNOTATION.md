# Variant annotation

## GRCh37

### SNVs and indels
Basic annotation of merged `vcf` files from the individual variants callers is carried out in two steps. First, the combined `vcf` is annotated with information from [RepeatMasker](http://www.repeatmasker.org/) and the [ENCODE consortium](http://rohsdb.cmb.usc.edu/GBshape/ENCODE/index.html). These files are retrieved from the [UCSC genome browser](https://genome.ucsc.edu) and parsed as such:

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

Subsequently, [vcf2maf](https://github.com/mskcc/vcf2maf) is used to annotate functional effects of mutations as well as other metadata using [VEP](https://www.ensembl.org/vep). The `--custom-enst` argument to vcf2maf takes a list of preferred gene transcript isoforms which to map mutations onto. We supply a consensus list of [`isoform_overrides_at_mskcc` and `isoform_overrides_uniprot`](https://github.com/mskcc/vcf2maf/tree/master/data), generated as such:
``` r
t1 = readr::read_tsv('isoform_overrides_at_mskcc')
t2 = readr::read_tsv('isoform_overrides_uniprot')
t2 %>%
    dplyr::filter(gene_name %nin% t1$gene_name) %>%
    dplyr::bind_rows(., t1) %>%
    readr::write_tsv('isoforms')
```