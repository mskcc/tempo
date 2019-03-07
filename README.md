## Vaporware pipeline

Pipeline for processing WES/WGS FASTQs into analysis-ready BAMs using GATK4 best practices, built using Nextflow domain specific language. 

See Nextflow documention here: 
https://www.nextflow.io/

Inspiration from `Sarek` at the NBIS and SciLifeLab in Sweden:
https://github.com/SciLifeLab/Sarek

Vaporware containes several pipeline scripts (alignment.nf, somatic.nf). Detalied explanation for each step is provided bellow.

### Setup instructions

Clone the repository as `vaporware`:

```
git clone https://github.com/mskcc/vaporware.git vaporware
cd vaporware
```

Install `Nextflow` within this subdirectory if the Nextflow executable `nextflow` isn't in your PATH:

```
curl -s https://get.nextflow.io | bash 
```

This installs Nextflow within the `vaporware` (the current directory). 

Optionally, move the `nextflow` file to a directory accessible by your $PATH variable in order to avoid remembering and typing the full path to `nextflow` each time.

#### Executing the scripts

Command for running nextflow

```
nextflow run <nextflow_script> --sample <samples_tsv>
```

nextflow_script - alignment.nf or somatic.nf
samples_tsv - samples - description how to create .tsv input files is provided below in steps section for each component of the pipeline

#### For submitting via LSF on juno

**NOTE:** In order to run successfully on `juno`, you must set `NXF_SINGULARITY_CACHEDIR`. This is not a random subdirectory. 

In order to get this working, use e.g. 
`export NXF_SINGULARITY_CACHEDIR=$HOME/.singularity`

See here for details: https://www.nextflow.io/docs/latest/config.html

```
NXF_SINGULARITY_CACHEDIR:
Directory where remote Singularity images are stored. When using a computing cluster it must be a shared folder accessible from all computing nodes.
```

* Do the following for LSF on juno:

```
nextflow run alignment.nf --sample test_inputs/lsf/test_samples.tsv -profile lsf_juno
```

#### For submitting via AWS Batch

In order to run pipeline on `AWS Batch`, you first must create `AWS Batch Compute Environment` and `AWS Batch Job Queue` as described here https://docs.aws.amazon.com/batch/latest/userguide/compute_environments.html, https://docs.aws.amazon.com/batch/latest/userguide/create-job-queue.html.
For user-specified Ami ID you should use `ami-077e66e85f2156f67`.
You must also create `S3 bucket` which will be used as working directory and where your outputs will be stored, as described here https://docs.aws.amazon.com/AmazonS3/latest/user-guide/create-bucket.html. 

**NOTE:** Your AWS Batch Instance role needs to have permissions to READ and WRITE to the S3 bucket.

Create `awsbatch.config` file in `conf` directory using `conf/awsbatch.config.template`

Replace:
`<AWS_REGION>` with aws region where you built your compute environment
`<STORAGE_ENCRYPTION>` storage encryption for your bucket (default: 'AES256')
`<AWS_BATCH_QUEUE_ARN>` ARN of your AWS Batch Job Queue
`<AWS_S3_WORKDIR>` S3 bucket used as working directory

* Do the following for AWS Batch:

```
nextflow run alignment.nf --sample test_inputs/aws/test_samples.tsv -profile awsbatch
```

#### Local, Docker, and Singularity

**NOTE:** To being able to run locally you need to provide reference files from `conf/references.config` and create `samples.tsv` as described in `Bioinformatic Components for the Alignment Script` section below.

The default parameters are for local use WITHOUT containers

* For Docker use, do the following:

```
nextflow run alignment.nf --sample samples.tsv -profile docker
```

* For Singularity use, do the following:

```
nextflow run alignment.nf --sample samples.tsv -profile singularity
```

## Components

### Bioinformatic Components for the Alignment Script

For the script `alignment.nf`, the pipeline does alignment with `bwa mem`, converts the SAM to a sorted BAM with `samtools`, and does uses `GATK4` to mark duplicates and do base recalibration. 

Input File columns:
`"patientId gender  status  sample  lane  fastq1  fastq2"`

Outputs:
They are found in `${params.outDir}/VariantCalling/<tool_name>`

Variables used in pipeline:

`patientId`: Patient Id

`gender`: XX or XY 

`status`: 0 - Normal 1 - Tumor

`sample`: Sample Id

`lane`: if sample is multiplexed

`fastq1`: Path to first pair of fastq

`fastq2`: Path to second pair of fastq

Execution on lsf:
```
nextflow run alignment.nf --sample test_inputs/lsf/test_alignment.tsv -profile lsf_juno --outDir $PWD
```

Execution on aws:
```
nextflow run alignment.nf --sample test_inputs/aws/test_alignment.tsv -profile awsbatch
```

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

### Bioinformatic Components for the Variant Calling Script

For the script `somatic.nf`, the pipeline performs the following:

```
delly call -> delly filter
snppileup -> facets
manta -> strelka
mutect2
```

NOTE: `manta -> strelka` and `mutect2` were lifted and modified from Sarek's implementation at https://github.com/mskcc/Sarek

Input File columns:
`"idTumor   idNormal    bamTumor    bamNormal   baiTumor    baiNormal"`

Outputs:
They are found in `${params.outDir}/VariantCalling/<tool_name>`

Variables used in pipeline:

`genomeFile`: reference fasta

`idTumor`: tumor sample name 

`idNormal`: normal sample name

`bamTumor`: tumor bam

`bamNormal`: normal bam


Execution on lsf:
```
nextflow run somatic.nf --sample test_inputs/lsf/test_somatic.tsv -profile lsf_juno --outDir $PWD
```

Execution on aws:
```
nextflow run somatic.nf --sample test_inputs/aws/test_somatic.tsv -profile awsbatch
```

You can also specifiy which specific tool(s) you want to run with the `--tools <toolname(s), comma-delimited>` flag.

If `--tools` is not specified, `somatic.nf` runs all tools, currently `delly`,`facets`,`manta`,`strelka2`, and `mutect2`.

Tool name `delly` will run processes `dellyCall`,`makeSamplesFile`, and `dellyFilter`.

Tool name `manta` runs ONLY process `runManta`; `manta` and `strelka2` needed to run process `runStrelka`.

Tool name `facets` runs process `doSNPPileup` and `doFacets`.

Tool name `mutect2` runs process `runMutect2`.

Example run of only `delly`, `manta, `strelka2` on lsf:
```
nextflow run somatic.nf --sample test_inputs/lsf/test_somatic.tsv -profile lsf_juno --outDir $PWD --tools delly, manta, strelka2
```

#### `delly` -- SV Caller 

https://github.com/dellytools/delly

For now we assume that we want all five variants of SV performed for `delly call`, thus the for-loop in the Nextflow process `dellyCall`.
```
sv_variants=("DUP" "BND" "DEL" "INS" "INV")
for sv_variant in "\${sv_variants[@]}";
do
  outfile="${idTumor}_${idNormal}_\${sv_variant}.bcf"
  delly call \
    -t "\${sv_variant}" \
    -o "\${outfile}" \
    -g ${genomeFile} \
    ${bamTumor} \
    ${bamNormal}
done
```

The next step, `delly filter`, requires a sample input file that contains a tumor ID and a normal ID that's considered "control"; this is done in the process `makeSampleFile`:

```
echo "${idTumor}\ttumor\n${idNormal}\tcontrol" > samples.tsv
```

This is fed to process `dellyFilter`; we also need a for-loop here for ease:

```
sv_variants=("DUP" "BND" "DEL" "INS" "INV")
for sv_variant in "\${sv_variants[@]}";
do
  delly_call_file="${idTumor}_${idNormal}_\${sv_variant}.bcf"
  outfile="${idTumor}_${idNormal}_\${sv_variant}.filter.bcf"
  delly filter \
    -f somatic \
    -o "\${outfile}" \
    -s "samples.tsv" \
    "\${delly_call_file}"
done
```

#### `MuTect2` -- SV Caller

https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_hellbender_tools_walkers_mutect_Mutect2.php

Currently our implementation of `mutect2` has hard-coded `-Xmx` parameters; should be adjusted later.

Also needs to be adjusted to work with scattered intervals.

```
gatk --java-options "-Xmx8g" \
  Mutect2 \
  -R ${genomeFile}\
  -I ${bamTumor}  -tumor ${idTumor} \
  -I ${bamNormal} -normal ${idNormal} \
  -O "${idTumor}_vs_${idNormal}_somatic.vcf"
```

#### `Manta` to `Strelka2` -- small variant caller

```
configManta.py \
--normalBam ${bamNormal} \
--tumorBam ${bamTumor} \
--reference ${genomeFile} \
--runDir Manta

python Manta/runWorkflow.py -m local -j ${task.cpus}

mv Manta/results/variants/candidateSmallIndels.vcf.gz \
  Manta_${idTumor}_vs_${idNormal}.candidateSmallIndels.vcf.gz
mv Manta/results/variants/candidateSmallIndels.vcf.gz.tbi \
  Manta_${idTumor}_vs_${idNormal}.candidateSmallIndels.vcf.gz.tbi
mv Manta/results/variants/candidateSV.vcf.gz \
  Manta_${idTumor}_vs_${idNormal}.candidateSV.vcf.gz
mv Manta/results/variants/candidateSV.vcf.gz.tbi \
  Manta_${idTumor}_vs_${idNormal}.candidateSV.vcf.gz.tbi
mv Manta/results/variants/diploidSV.vcf.gz \
  Manta_${idTumor}_vs_${idNormal}.diploidSV.vcf.gz
mv Manta/results/variants/diploidSV.vcf.gz.tbi \
  Manta_${idTumor}_vs_${idNormal}.diploidSV.vcf.gz.tbi
mv Manta/results/variants/somaticSV.vcf.gz \
  Manta_${idTumor}_vs_${idNormal}.somaticSV.vcf.gz
mv Manta/results/variants/somaticSV.vcf.gz.tbi \
  Manta_${idTumor}_vs_${idNormal}.somaticSV.vcf.gz.tbi
```

`mantaCSI` is `Manta_${idTumor}_vs_${idNormal}.candidateSmallIndels.vcf.gz`, made in `manta` step and passed to `strelka2`

```
configureStrelkaSomaticWorkflow.py \
--tumor ${bamTumor} \
--normal ${bamNormal} \
--referenceFasta ${genomeFile} \
--indelCandidates ${mantaCSI} \
--runDir Strelka

python Strelka/runWorkflow.py -m local -j ${task.cpus}

mv Strelka/results/variants/somatic.indels.vcf.gz \
  Strelka_${idTumor}_vs_${idNormal}_somatic_indels.vcf.gz
mv Strelka/results/variants/somatic.indels.vcf.gz.tbi \
  Strelka_${idTumor}_vs_${idNormal}_somatic_indels.vcf.gz.tbi
mv Strelka/results/variants/somatic.snvs.vcf.gz \
  Strelka_${idTumor}_vs_${idNormal}_somatic_snvs.vcf.gz
mv Strelka/results/variants/somatic.snvs.vcf.gz.tbi \
  Strelka_${idTumor}_vs_${idNormal}_somatic_snvs.vcf.gz.tbi
```

#### `snp-pileup` to `doFacets` -- CNV caller

`doFacets.R` requires a counts file, so `snp-pileup` is ran first in process `doSNPPileup`.

```
output_filename = idTumor + "_" + idNormal + ".snppileup.dat.gz"
snp-pileup -A -P 50 --gzip "${facetsVcf}" "${output_filename}" "${bamTumor}" "${bamNormal}"
```

Process `doFacets` needs parameters generalized, but for now it is shown below:

```
snp_pileup_prefix = idTumor + "_" + idNormal
counts_file = "${snp_pileup_prefix}.snppileup.dat.gz"
genome_value = "hg19"
TAG = "${snp_pileup_prefix}"
directory = "."

/usr/bin/facets-suite/doFacets.R \
--cval 100 \
--snp_nbhd 250 \
--ndepth 35 \
--min_nhet 25 \
--purity_cval 500 \
--purity_snp_nbhd 250 \
--purity_ndepth 35 \
--purity_min_nhet 25 \
--genome "${genome_value}" \
--counts_file "${counts_file}" \
--TAG "${TAG}" \
--directory "${directory}" \
--R_lib latest \
--single_chrom F \
--ggplot2 T \
--seed 1000 \
--tumor_id "${idTumor}"
```
