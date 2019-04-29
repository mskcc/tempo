# Vaporware 

CMO WES/WGS pipeline using [Nextflow framework](https://github.com/nextflow-io/nextflow)

See Nextflow documention here: 
https://www.nextflow.io/

Inspiration from `Sarek` at the NBIS and SciLifeLab in Sweden:
https://github.com/SciLifeLab/Sarek

Vaporware contains several pipeline scripts (`make_bam_and_qc.nf`, `somatic.nf`, `germline.nf`). A detalied explanation for each step is provided below.

## Setup instructions

Clone the repository as `vaporware`:

```
git clone https://github.com/mskcc/vaporware.git vaporware
cd vaporware
git submodule update --recursive --remote    ## pull all remote submodules
```

Install `Nextflow` within this subdirectory if the Nextflow executable `nextflow` isn't in your PATH:

```
curl -s https://get.nextflow.io | bash 
```

This installs Nextflow within the subdirectory `vaporware`.

Optionally, move the `nextflow` file to a directory accessible by your `$PATH` variable in order to avoid remembering and typing the full path to `nextflow` each time.

## Executing the scripts

Command for running nextflow

```
nextflow run <nextflow_script> --sample <samples_tsv>
```

* `<nextflow_script>` - e.g. `make_bam_and_qc.nf` or `somatic.nf`

* `<samples_tsv>` - sample file

### For submitting via LSF on juno

**NOTE:** In order to run successfully on `juno`, you must set `NXF_SINGULARITY_CACHEDIR`. Use e.g. 
`export NXF_SINGULARITY_CACHEDIR=$HOME/.singularity`

See here for details: https://www.nextflow.io/docs/latest/config.html

We recommend users use `singularity/3.1.1`, i.e. 

```
$ module load  singularity/3.1.1
$ which singularity
/opt/local/singularity/3.1.1/bin/singularity
```

```
NXF_SINGULARITY_CACHEDIR:
Directory where remote Singularity images are stored. 
When using a computing cluster it must be a shared folder accessible from all computing nodes.
```

* Do the following for LSF on juno:

```
nextflow run make_bam_and_qc.nf --mapping test_inputs/lsf/test_make_bam_and_qc.tsv --pairing test_inputs/lsf/test_make_bam_and_qc_pairing.tsv -profile juno
```

### For submitting via AWS Batch

In order to run pipeline on `AWS Batch`, you first must create your `Compute Environment` and `Configuration File` as described [here](aws_cf_scripts/README.md).
 
* When you build your compute environment and create configuration file do the following:

```
nextflow run make_bam_and_qc.nf --mapping test_inputs/aws/test_make_bam_and_qc.tsv --pairing test_inputs/aws/test_make_bam_and_qc_pairing.tsv -profile awsbatch
```

### Local, Docker, and Singularity

**NOTE:** To being able to run locally you need to provide reference files from `conf/references.config` and create `samples.tsv` as described in `Bioinformatic Components for the Make Bam and QC Script` section below.

The default parameters are for local use WITHOUT containers

* For Docker use, do the following:

```
nextflow run make_bam_and_qc.nf --mapping mapping.tsv --pairing pairing.tsv -profile docker
```

* For Singularity use, do the following:

```
nextflow run make_bam_and_qc.nf --mapping mapping.tsv --pairing pairing.tsv -profile singularity
```

## Components

### Bioinformatic Components for the Make Bam and QC Script

For the script `make_bam_and_qc.nf`, the pipeline does alignment with `bwa mem`, converts the SAM to a sorted BAM with `samtools`, and does uses `GATK4` to mark duplicates and do base recalibration. 

Input File columns:
`"patientId gender  status  sample  lane  fastq1  fastq2"`

Outputs:
They are found in `${params.outDir}/VariantCalling/<tool_name>`

Variables used in pipeline:

### Mapping file:

`sample`: Sample Id

`lane`: if sample is multiplexed

`assay`: exome or genome

`target`: target file

`fastq1`: Path to first pair of fastq

`fastq2`: Path to second pair of fastq

### Pairing file:

`normalId`: normal sampleId

`tumorId` : tumor sampleId

Execution on lsf:
```
nextflow run make_bam_and_qc.nf --mapping test_inputs/lsf/test_make_bam_and_qc.tsv --pairing test_inputs/lsf/test_make_bam_and_qc_pairing.tsv -profile juno --outDir $PWD
```

Execution on aws:
```
nextflow run make_bam_and_qc.nf --mapping test_inputs/aws/test_make_bam_and_qc.tsv --pairing test_inputs/aws/test_make_bam_and_qc_pairing.tsv -profile awsbatch
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

 ### Outputs:

 `make_bam_and_qc.nf` outputs `make_bam_output.tsv` which can be used as input for `somatic.nf`. You can change the name of the output file using `--outname filename.tsv` parametar.

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
nextflow run somatic.nf --sample test_inputs/lsf/test_somatic.tsv -profile juno --outDir $PWD
```

Execution on aws:
```
nextflow run somatic.nf --sample test_inputs/aws/test_somatic.tsv -profile awsbatch
```

You can also specifiy which specific tool(s) you want to run with the `--tools <toolname(s), comma-delimited>` flag.

If `--tools` is not specified, `somatic.nf` runs all tools, currently `delly`,`facets`,`manta`,`strelka2`, `mutect2`, and `msisensor`.

Tool name `delly` will run processes `dellyCall`,`makeSamplesFile`, and `dellyFilter`.

Tool name `manta` runs ONLY process `runManta`; `manta` and `strelka2` needed to run process `runStrelka`.

Tool name `facets` runs process `doSNPPileup` and `doFacets`.

Tool name `mutect2` runs process `runMutect2`.

Tool name `msisensor` runs process `runMsiSensor`.

Example run of only `delly`, `manta, `strelka2` on lsf:
```
nextflow run somatic.nf --sample test_inputs/lsf/test_somatic.tsv -profile juno --outDir $PWD --tools delly, manta, strelka2
```

Input File columns:
`"sequenceType  idTumor   idNormal    bamTumor    bamNormal   baiTumor    baiNormal"`

Outputs:
They are found in `${params.outDir}/VariantCalling/${idTumor}_${idNormal}/<tool_name>`

Variables used in pipeline:

`genomeFile`: reference fasta

`sequenceType`: Either `exome` or `genome`; currently un-used

`idTumor`: tumor sample name 

`idNormal`: normal sample name

`bamTumor`: tumor bam

`bamNormal`: normal bam

#### `delly` -- SV Caller 

https://github.com/dellytools/delly

For now we assume that we want all five variants of SV performed for `delly call`; `sv_variant` is defined before the process is initiated in a list named `sv_variants`. 

We use `nextflow` to iterate through that list, submitting the processes in parallel.

```
outfile="${idTumor}_${idNormal}_\${sv_variant}.bcf"
delly call \
  -t "\${sv_variant}" \
  -o "\${outfile}" \
  -g ${genomeFile} \
  ${bamTumor} \
  ${bamNormal}
```

The next step, `delly filter`, requires a sample input file that contains a tumor ID and a normal ID that's considered "control"; this is done in the process `makeSampleFile`:

```
echo "${idTumor}\ttumor\n${idNormal}\tcontrol" > samples.tsv
```

These `sv_variant` processes are then fed to another process `dellyFilter`.

```
delly_call_file="${idTumor}_${idNormal}_\${sv_variant}.bcf"
outfile="${idTumor}_${idNormal}_\${sv_variant}.filter.bcf"
delly filter \
  -f somatic \
  -o "\${outfile}" \
  -s "samples.tsv" \
  "\${delly_call_file}"
```

#### `MuTect2` -- SV Caller

https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_hellbender_tools_walkers_mutect_Mutect2.php

Running `mutect2` consists of three processes: `CreateIntervalBeds`, `runMutect2`, `indexVCF`, `runMutect2Filter`, and `combineMutect2VCF`.

The output of `combineMutect2VCF` is later then sent to process `combineChannel` to be merged with `strelka2` output into one MAF file.

##### SplitIntervals

Before running `mutect2`, we use `SplitIntervals` (in process `CreateIntervalBeds`) against interval file `human.b37.genome.bed`. Currently we only split into a `scatter-count` of 10; we can adjust in the future.

```
  gatk SplitIntervals \
    -R ${genomeFile} \
    -L ${intervals} \
    --scatter-count 10 \
    --subdivision-mode BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW \
    -O interval_beds
```

##### runMutect2

Currently our implementation of process `runMutect2` has hard-coded `-Xmx` parameters; should be adjusted later.

NOTE: A somatic.nf run of `mutect2` will create a certain number of processes for each of the subsequent steps based on how many splits are performed above (currently set at 10), up until `combineMutect2VCF`.

```
gatk --java-options "-Xmx8g" \
  Mutect2 \
  -R ${genomeFile}\
  -I ${bamTumor}  -tumor ${idTumor} \
  -I ${bamNormal} -normal ${idNormal} \
  -O "${idTumor}_vs_${idNormal}_somatic.vcf"
```

##### Index, Filter, and Combine

The outputs from `runMutect2` - `mutect2Vcf` is indexed with `tabix` (from process `indexVCF`):

```
  tabix -p vcf ${mutect2Vcf}
```

This creates a `.tbi` file stored in `mutect2VcfIndex`. They are then fed into `FilterMutectCalls` (from process `runMutect2Filter`):

```
  # Xmx hard-coded for now due to lsf bug
  gatk --java-options "-Xmx8g" \
    FilterMutectCalls \
    --variant "${mutect2Vcf}" \
    --output "${outfile}"
```

All generated, filtered `mutect2` `vcf` files are then submitted to process `combineMutect2VCF`. This uses `bcftools concat | bcftools sort` to create the final `mutect2` output VCF.

```
  outfile="${idTumor}_${idNormal}.mutect2.filtered.combined.vcf.gz"

  bcftools concat ${mutect2Vcfs} | bcftools sort --output-type z --output-file ${outfile}
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

#### combineChannel; bcftools filter, norm, and merge

The `combineChannel` process just takes the VCF files from the `combineMutect2VCF` and `runStrelka` processes into one Channel, `vcfOutputSet`. This allows parallel job submission for each VCF.

Process `runBCFToolsFilterNorm` performs three commands to the `${vcf}` file: `tabix` to index it; `bcftools filter` to filter the file and output a zipped `vcf`; and `bcftools norm` to perform a normalization step.

```
  tabix -p vcf ${vcf}

  bcftools filter \
    -r 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,MT,X,Y \
    --output-type z \
    "${vcf}" | \
  bcftools norm \
    --fasta-ref ${genomeFile} \
    --output-type z \
    --output "${outfile}"
```

#### bcftools merge

The process `runBCFToolsMerge` loads `vcfFilterNormOutput` as input.

Each `${vcf}` in `vcfFilterNormOutput` is first indexed with `tabix` and then fed into `bcftools merge`. The output is an uncompressed VCF file so that it can be processed by `vcf2maf`.

Output is sent to channel `vcfMergedOutput`.

```
  for f in *.vcf.gz
  do
    tabix -p vcf \$f
  done

  bcftools merge \
    --force-samples \
    --merge none \
    --output-type v \
    --output "${idTumor}_${idNormal}.mutect2.strelka2.filtered.norm.merge.vcf" \
    *.vcf.gz
```

#### vcf2maf

Process `runVCF2MAF` takes the VCF file in channel `vcfMergedOutput` and runs `vcf2maf` to make one MAF file.

The hard-coded paths in the below script are relative to the container used by the process.

```
  perl /opt/vcf2maf.pl \
    --input-vcf ${vcfMerged} \
    --tumor-id ${idTumor} \
    --normal-id ${idNormal} \
    --vep-path /opt/vep/src/ensembl-vep \
    --vep-data ${vepCache} \
    --filter-vcf ${vcf2mafFilterVcf} \
    --output-maf ${outfile} \
    --ref-fasta ${genomeFile}
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

#### `msisensor` -- MicroSatellite Instability detection

https://github.com/ding-lab/msisensor

msiSensorList is defined in `references.config`

```
output_prefix = "${idTumor}_${idNormal}"

msisensor msi -d "${msiSensorList}" -t "${bamTumor}" -n "${bamNormal}" -o "${output_prefix}"
```

