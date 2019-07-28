# Reference files

## GRCh37

### Genome
[GATK bundle](https://software.broadinstitute.org/gatk/download/bundle), also available [here](https://console.cloud.google.com/storage/browser/gatk-legacy-bundles/b37).

### SNV and indel calling
For genomes, use [`bed` file of "callable"](https://github.mskcc/vaporware/blob/master/docs/INTERVALS.md#genome) regions from GATK's bundle.

For exomes, use `bed` file corresponding to the platform used for target capture, see [documentation on intervals](https://github.mskcc/vaporware/blob/master/docs/INTERVALS.md#exome-capture-platform).


### Structural variant calling
Delly provides and takes as an argument a [file of regions](https://github.com/dellytools/delly/tree/master/excludeTemplates) to _exclude_ from variant calling. This excludes telomeres and centromeres from auto- and allosomes as well as any other contig.

For Manta, substract these regions from a bed file of the whole genome to generate a list of regions to _include_. First clean up the file provided by Delly, since it is not in `bed` format:
``` shell
grep -Ev "chr|MT|GL00|NC|hs37d5" human.hg19.excl.tsv > human.hg19.excl.clean.bed
bedtools subtract -a b37.bed -b human.hg19.excl.clean.bed > b37.minusDellyExclude.bed
```

## GRCh38
[GATK bundle](https://software.broadinstitute.org/gatk/download/bundle), also available [here](https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38).

**Note:** Support for hg38 is currently somewhat limited. Please raise an issue at https://github.com/mskcc/vaporware/issues if you would like to process data with GRch38. However, please note that hg38 reference files are easily avilable from UCSC. 

