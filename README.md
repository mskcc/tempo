# vaporware
pipeline development repo

## Tools, Docker images, and Files used in CWL Pipeline Example
### Based on [GATK Best Practices](https://software.broadinstitute.org/gatk/best-practices/workflow?id=11165)

**FASTQ gz test file**

Grab `P1_R1.fastq.gz` and `P1_R2.fastq.gz` in the Roslin github repo `https://github.com/mskcc/roslin-variant/tree/master/setup/examples/data/fastq/`

**Reference Files**

*Reference files located in luna or selene*

1. FASTA GRCh37 Reference Genome and indexes (.amb, .ann, .fai, etc) from `/ifs/depot/pi/resources/genomes/GRCh37/fasta/b37.*`

2. Sequence group interval data `/ifs/depot/pi/resources/targets/IMPACT468_b37/targets_list/picard_targets.interval_list`

3. VCF files and indexes (don't forget to grab the .idx)
```
/ifs/depot/pi/resources/genomes/GRCh37/indels_1000g/Mills_and_1000G_gold_standard.indels.b37.vcf*
/ifs/depot/pi/resources/genomes/GRCh37/hapmap33/hapmap_3.3.b37.vcf*
/ifs/depot/pi/resources/genomes/GRCh37/dbsnp/138/dbsnp_138.b37.vcf*
```

*On AWS*

Example data in `s3:testbucketcmo/example_data/fastq/`.
Reference files in `s3:testbucketcmo/ref/grch37/`.

**Alignment**

1. bwa-mem (quay.io/collaboratory/dockstore-tool-bwa-mem:1.0). bwa-mem also needs a `read_group` parameter, or else downstream steps don't work; you can use `"@RG\\tID:Seq01p\\tSM:Seq01\\tPL:ILLUMINA\\tPI:330"`

2. samtools view (quay.io/cancercollaboratory/dockstore-tool-samtools-view:1.0); this converts bwa-mem sam output to bam

3. samtools sort (quay.io/cancercollaboratory/dockstore-tool-samtools-sort:1.0)

**Mark Duplicates**
1. gatk MarkDuplicates (broadinstitute/gatk:latest)

**BQSR**
1. gatk BuildBamIndex (broadinstitute/gatk:latest)

2. gatk BaseRecalibrator (broadinstitute/gatk:latest)

3. gatk ApplyBQSR (broadinstitute/gatk:latest)
