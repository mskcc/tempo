# Genomic intervals

## GRCH37

### Genome
Create a `bed` of the auto- and allosomes in b37 (refer to `/juno/work/taylorlab/cmopipeline/mskcc-igenomes/grch37/genome/b37.bed`) by downloading chromosome sizes and parsing the output into proper format:
``` shell
wget https://raw.githubusercontent.com/igvteam/igv/master/genomes/sizes/b37.chrom.sizes
paste <(cut -f1 b37.chrom.sizes | head -24) \
      <(seq 24 | awk '{print "0"}') \
      <(cut -f2 b37.chrom.sizes | head -24) \
      > b37.bed

cat b37.bed
1       0       249250621
2       0       243199373
3       0       198022430
4       0       191154276
5       0       180915260
6       0       171115067
7       0       159138663
8       0       146364022
9       0       141213431
10      0       135534747
11      0       135006516
12      0       133851895
13      0       115169878
14      0       107349540
15      0       102531392
16      0       90354753
17      0       81195210
18      0       78077248
19      0       59128983
20      0       63025520
21      0       48129895
22      0       51304566
X       0       155270560
Y       0       59373566
```

Generate the `bed` file of "callable" regions as such:
``` shell
gatk IntervalListToBed --INPUT b37_wgs_calling_regions.v1.interval_list --OUTPUT b37_wgs_calling_regions.v1.bed
```

### Exome Capture Platform
Currently supporting:
- __AgilentExon_51MB__: SureSelectXT Human All Exon V4 from Agilent
- __IDT_Exome__: xGen Exome Research Panel v1.0 from IDT

Add 5 bp to each end of exons to make sure splice site mutations can be called:
``` shell
bedtools slop \
    -g b37.chrom.sizes \
    -i targets.bed \
    -r 5 \
    -l 5 \
    > targets.plus5bp.bed
```

## GRCh38
[GATK bundle](https://software.broadinstitute.org/gatk/download/bundle), also available [here](https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38).

**Note:** Support for hg38 is currently somewhat limited. Please raise an issue at https://github.com/mskcc/vaporware/issues if you would like to process data with GRch38. However, please note that hg38 reference files are easily avilable from UCSC. 

