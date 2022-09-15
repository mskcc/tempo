# Reference Files

This and associated pages in this section provide details on the provenance and generation of all reference files used in `pipeline.nf`. Usage of these files is defined in the [references configuration file](https://github.com/mskcc/tempo/blob/master/conf/references.config).

::: tip Note
All reference files described herein are in assembly GRCh37/hg19 of the human genome.
:::

## Genome Assembly

Part of the [GATK bundle](https://software.broadinstitute.org/gatk/download/bundle), also available [here](https://console.cloud.google.com/storage/browser/gatk-legacy-bundles/b37). Tempo uses the **human_g1k_v37_decoy** assembly of the genome.

## Genomic Intervals

BED files that specify the regions of the genome to consider for variant calling are specified in the [input files](running-the-pipeline.md#input-files).

### Exome Capture Platforms
For exomes, use BED file corresponding to the platform used for target capture. Currently, Juno reference files are configured to support:
- __AgilentExon_51MB__: SureSelectXT Human All Exon V4 from Agilent.
- __IDT_Exome__: xGen Exome Research Panel v1.0 from IDT.
- __IDT_Exome_v2__: xGen Exome Research Panel v2.0 from IDT.

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

## Custom target files

Required target files have already been built for Agilent and IDT exome baits and can be used readily with Tempo. Read this section if you have a different target design in mind. 

### Required files

* `targets.bed`: The targets file can be obtained from the provider of the baitset. If you only have the file in gzipped format, please run `zcat <targets>.bed.gz > <targets>.bed`. In practice, we alter this file to have 5 bp padding on each region in both directions. Padding can be added using the `slop` subcommand of `bedtools`. 
* `targets.bed.gz`: You also need a gzipped copy of the `targets.bed` file. Create with `bgzip`
* `targets.bed.gz.tbi`: The above file will need to be indexed with `tabix`.
* `targets.interval_list`: The bed file can be used as input to create the interval list. Create using the `BedToIntervalList` tools from gatk.
* `baits.interval_list`: The baits file can also be obtained from the provider of the baitset, typically as a bed file. Create the interval list using the `BedToIntervalList` tool from gatk. If a bait bed file is not provided, you can copy or link to `targets.interval_list` instead.
* `coding.bed`: This file will be used to calculate TMB. The known coding regions of the reference should first be downloaded (ex: [EnsGene for hg19](https://genome.ucsc.edu/cgi-bin/hgTables?hgsid=1138949195_54NfeOmJerLbPvAdqAA6vaGWonRr&clade=mammal&org=Human&db=hg19&hgta_group=genes&hgta_track=ensGene&hgta_table=0&hgta_regionType=genome&position=chrX%3A15%2C578%2C261-15%2C621%2C068&hgta_outputType=bed&hgta_outFileName=)) and then be intersected with the targets file using `bedtools`. The result should subsequently be sorted and merged with `bedtools`. 

### Input to Tempo

You should designate a folder just for the required file. This folder will be input as a parameter called `targets_base`, and all of your targets should be placed there. The folder structure should looks as follows:
``` shell
<targets_base folder>
├── <target name 1>
├── <target name 2>
└── <target name 3>
```
You can have as many targets as you like. Under each folder are the 6 target files previously described. For example:
``` shell
<targets_base folder>/agilent/
├── baits.interval_list -> ../../baits/AgilentExon_51MB_b37_v3_baits.interval_list
├── coding.bed -> ../../coding_regions/AgilentExon_51MB_b37_v3_baits.coding.sorted.merged.bed
├── targets.bed -> ../../targets/AgilentExon_51MB_b37_v3_targets.plus5bp.bed
├── targets.bed.gz -> ../../targets/AgilentExon_51MB_b37_v3_targets.plus5bp.bed.gz
├── targets.bed.gz.tbi -> ../../targets/AgilentExon_51MB_b37_v3_targets.plus5bp.bed.gz.tbi
└── targets.interval_list -> ../../targets/AgilentExon_51MB_b37_v3_targets.interval_list
```
In this case, the files are symbolically linked to the original. Whether soft links or hard links are used, the files in this folder should strictly match the names `coding.bed`, `baits.interval_list`, `targets.interval_list`, `targets.bed`, `targets.bed.gz`, `targets.bed.gz.tbi`. 

When running Tempo, use the parameter `--targets_base <targets_base folder>` so that Nextflow will know where to find your target files. 
When running with `--assayType genome`, only the `<targets_base>/wgs` target folder will be available. Conversely, the `<targets_base>/wgs` target folder will not be available when `--assayType genome` is not set.

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

For BRASS, Tempo is using a pre-built reference packages linked to the [dockstore registry](ftp://ftp.sanger.ac.uk/pub/cancer/dockstore/human/). To download references for GRCh37:
``` shell
wget ftp://ftp.sanger.ac.uk/pub/cancer/dockstore/human/VAGrENT_ref_GRCh37d5_ensembl_75.tar.gz
tar -xzvf VAGrENT_ref_GRCh37d5_ensembl_75.tar.gz
wget ftp://ftp.sanger.ac.uk/pub/cancer/dockstore/human/CNV_SV_ref_GRCh37d5_brass6+.tar.gz
tar -xzvf CNV_SV_ref_GRCh37d5_brass6+.tar.gz
# grab CNV_SV_ref/brass and VAGrENT_ref_GRCh37d5_ensembl_75/vagrent
```
GRCh38 is also available from the same ftp server. To build a new reference, follow the instructions on the [BRASS wiki](https://github.com/cancerit/BRASS/wiki).

## Structural Variant Annotation

The bed and bedpe files used for the flags `pcawg_blacklist_bed`,`pcawg_blacklist_bedpe`, `pcawg_blacklist_fb_bedpe` and `pcawg_blacklist_te_bedpe` are sourced from the [SV merging tool used in the PCAWG paper](https://bitbucket.org/weischenfeldt/pcawg_sv_merge/src/docker/data/blacklist_files/). They can be downloaded as follows: 
``` shell
wget https://api.bitbucket.org/2.0/repositories/weischenfeldt/pcawg_sv_merge/src/docker/data/blacklist_files/pcawg6_blacklist.slop.bed.gz
wget https://api.bitbucket.org/2.0/repositories/weischenfeldt/pcawg_sv_merge/src/docker/data/blacklist_files/pcawg6_blacklist.slop.bedpe.gz
wget https://api.bitbucket.org/2.0/repositories/weischenfeldt/pcawg_sv_merge/src/docker/data/blacklist_files/pcawg6_blacklist_foldback_artefacts.slop.bedpe.gz
wget https://api.bitbucket.org/2.0/repositories/weischenfeldt/pcawg_sv_merge/src/docker/data/blacklist_files/pcawg6_blacklist_TE_pseudogene.bedpe.gz
```

Structural variants are also filtered with RepeatMasker and Mappability blacklists, whose preparation are described in the [above section](#repeatMasker-and-mappability-blacklist).
