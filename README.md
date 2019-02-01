# vaporware

## Nextflow framework

Basic Nextflow pipline for processing WES/WGS FASTQs into analysis-ready BAMs using GATK4 best practices. The following pipeline should do:

* alignment via `bwa mem`
* `MarkDuplicates` to indentify duplicate reads for QC with GATK4
* Base recalibration via `gatk BaseRecalibrator` and `gatk ApplyBQSR`

See Nextflow documention here: 
https://www.nextflow.io/

Inspiration from `Sarek` at the NBIS and SciLifeLab in Sweden:
https://github.com/SciLifeLab/Sarek

### Usage instructions

Clone the repository as `vaporwareNextflow`:

```
git clone -b feature/nextflow https://github.com/mskcc/vaporwarevaporwareNextflow
cd vaporwareNextflow
```

(You need to use the branch `feature/nextflow` if cloning from master.)

Install `Nextflow` within this subdirectory if the Nextflow executable `nextflow` isn't globallly defined:

```
curl -s https://get.nextflow.io | bash 
```

This installs Nextflow within the `vaporwareNextflow` (the current directory)

###### Using the pipeline

```
nextflow run SciLifeLab/Sarek/main.nf --sample samples_germline.tsv
```

more details here...


### Bioinformatic Components

* `bwa mem` -- alignment

http://bio-bwa.sourceforge.net/bwa.shtml

```
bwa mem ref.fa reads.fq > aln-se.sam
bwa mem ref.fa read1.fq read2.fq > aln-pe.sam
```

* `gatk MarkDuplicates` --- Identifies duplicate reads 


https://software.broadinstitute.org/gatk/documentation/tooldocs/4.0.4.0/picard_sam_markduplicates_MarkDuplicates.php

```
gatk  MarkDuplicates \ 
    I=input.bam \ 
    O=marked_duplicates.bam \ 
    M=marked_dup_metrics.txt
```

* `gatk BaseRecalibrator` --- Detect systematic errors in base quality scores

https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_bqsr_BaseRecalibrator.php

```
 gatk BaseRecalibrator \
   -R reference.fasta \
   -I my_reads.bam \
   -knownSites latest_dbsnp.vcf \
   -o recal_data.table
 ```


* `gatk ApplyBQSR` --- Apply base quality score recalibration

https://software.broadinstitute.org/gatk/documentation/tooldocs/4.beta.2/org_broadinstitute_hellbender_tools_walkers_bqsr_ApplyBQSR.php

```
 gatk ApplyBQSR \
   -R reference.fasta \
   -I input.bam \
   -BQSR recalibration.table \
   -O output.bam
 
 ```