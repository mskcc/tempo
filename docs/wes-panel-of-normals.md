# Creating a Panel of Normals (PoN) for exomes

"Somatic" variants that occur in a panel of normal samples can be considered sequencing artifacts. We can generate a VCF file to filter against by calling variants in normal samples that look "clean", i.e. absent of tumor contamination. We use a similar variant calling strategy as for the somatic variant calling in tumor samples

For each normal sample call variants with `Strelka2` and `MuTect2`.

## Strelka2

Run `Manta` to seed indel calling, otherwise run as if the normal sample is an unmatched tumor sample. Parse output with `bcftools`, subsetting on variants supported by more than one alternate read.

```shell
$MANTA/configManta.py \
    --referenceFasta $REF \
    --runDir pon/manta/$NORMAL_NAME \
    --exome \
    --callRegions $TARGETS
    --bam $NORMAL_BAM

pon/manta/$NORMAL_NAME/runWorkflow.py --mode local

$STRELKA/configureStrelkaGermlineWorkflow.py \
    --ref $REF \
    --runDir mutations/pon/strelka2/$NORMAL_NAME \
    --exome \
    --callRegions $TARGETS \
    --indelCandidates pon/manta/$NORMAL_NAME/results/variants/candidateSmallIndels.vcf.gz \
    --bam $NORMAL_BAM

bcftools filter \
    --include 'FORMAT/AD[0:1]>1' \
    pon/strelka2/$NORMAL_NAME/results/variants/variants.vcf.gz | \
    bcftools norm \
    --fasta-ref $REF \
    -check-ref s \
    --multiallelics -both \
    --output-type z \
    --output pon/$NORMAL_NAME.strelka2.vcf.gz

tabix --preset vcf pon/$NORMAL_NAME.strelka2.vcf.gz
```

## MuTect2

`MuTect2` provides a variant calling mode for normal samples for this purpose. Process the output similarly to above. Fix some VCF header tags so that the files can be combined downstream. As opposed to the somatic variant calling in tumor samples, here retain any calls at multiallelic loci.

```shell
$GATK Mutect2 \
    --reference $REF \
    --intervals $TARGETS \
    --input $NORMAL_BAM \
     --tumor $NORMAL_NAME \
     --output pon/mutect2/$NORMAL_NAME.vcf.gz

bcftools filter \
    --include 'FORMAT/AD[0:1]>1' \
    pon/mutect2/$NORMAL_NAME.vcf.gz | \
    sed -e 's/ID=RU,Number=1/ID=RU,Number=A/' -e 's/ID=AD,Number=R/ID=AD,Number=./' |
    bcftools norm \
    --fasta-ref $REF \
    --check-ref s \
    --multiallelics -both \
    --output-type z \
    --output pon/$NORMAL_NAME.mutect2.vcf.gz

tabix --preset vcf pon/$NORMAL_NAME.mutect2.vcf.gz
```

Now, combine all individual VCFs from all normal samples. This requires a [`bcftools` plugin](https://samtools.github.io/bcftools/howtos/plugins.html).

```shell
bcftools merge \
    --merge none \
    --output-type z \
    --output pon.vcf.gz \
    pon/*vcf.gz

bcftools +fill-tags pon.vcf.gz \
    --output-type z \
    --output pon.annot.vcf.gz \
    -- --tags AC

tabix --preset vcf pon.annot.vcf.gz
```

Now, `pon.annot.vcf.gz` is ready to use to [annotate somatic variant calls from tumor samples](https://github.com/mskcc/vaporware/blob/135719e430b7e7338a1aff25831968b97267b343/pipeline.nf#L931).