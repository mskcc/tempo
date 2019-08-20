# gnomAD
<small>As of 2019-08-16, the latest release of gnomAD is 2.1.1 from March 6, 2019.</small>

Variant-level allele counts and frequencies from gnomAD's exome and genome cohorts are used to annotate somatic and germline SNVs/indels.

These files were retrieved and processed as such:

## Exomes

The exome VCF file contains a lot of information. We prune most of this to reduce file size and only keep relevant information. Only values from the non-cancer subset of the total population are used, since this excludes the normals from TCGA.

```shell
# Download files
wget https://storage.googleapis.com/gnomad-public/release/2.1.1/vcf/exomes/gnomad.exomes.r2.1.1.sites.vcf.bgz
wget https://storage.googleapis.com/gnomad-public/release/2.1.1/vcf/exomes/gnomad.exomes.r2.1.1.sites.vcf.bgz.tbi
mv gnomad.exomes.r2.1.1.sites.vcf.bgz gnomad.exomes.r2.1.1.sites.vcf.gz
mv gnomad.exomes.r2.1.1.sites.vcf.bgz.tbi gnomad.exomes.r2.1.1.sites.vcf.gz.tbi

# Parse INFO columns to retain only relevant ones
bcftools view --header-only gnomad.exomes.r2.1.1.sites.vcf.gz | \
    grep -E "non_cancer_AC|non_cancer_AF" | \
    grep -v -e "_male" -e "_female" \
    > gnomad.exomes.r2.1.1.sites.retained.info

paste <(grep -oP "(?<=ID\=)[A-Za-z_]+" gnomad.exomes.r2.1.1.sites.retained.info) \
      <(grep -oP "(?<=Description=\")[A-Za-z\(\),\-\_\ ]+" gnomad.exomes.r2.1.1.sites.retained.info) \
      > tmp && \
      mv tmp gnomad.exomes.r2.1.1.sites.retained.info

COLS=$(Rscript -e "out = paste0('^INFO/', paste(read.delim('gnomad.exomes.r2.1.1.sites.retained.info', header = F)[['V1']], collapse = ',INFO/')); cat(out)")

# Apply this, and retain filtered sites
bcftools annotate \
    --remove "$COLS" \
    --include 'FILTER~"PASS" | FILTER~"RF"' \
    --output-type z \
    --output tmp.vcf.gz \
    gnomad.exomes.r2.1.1.sites.vcf.gz

tabix --preset vcf tmp.vcf.gz

# Mark filtered sites
bcftools annotate \
    --annotations gnomad.exomes.r2.1.1.sites.non_cancer.vcf.gz \
    --include 'FILTER!="PASS"' \
    --mark-sites "+gnomAD_FILTER" \
    -k \
    --output-type z \
    --output gnomad.exomes.r2.1.1.sites.non_cancer.vcf.gz \
    tmp.vcf.gz

tabix --preset vcf gnomad.exomes.r2.1.1.sites.non_cancer.vcf.gz

# Clean up
rm gnomad.exomes.r2.1.1.sites.retained.info
rm tmp.vcf.gz tmp.vcf.gz.tbi
```

## Genomes

For genomes, there is no non-cancer subset.

```shell
# Download, one chromosome at the time
for chr in {1..22} X
do
    wget https://storage.googleapis.com/gnomad-public/release/2.1.1/vcf/genomes/gnomad.genomes.r2.1.1.sites.${chr}.vcf.bgz
    wget https://storage.googleapis.com/gnomad-public/release/2.1.1/vcf/genomes/gnomad.genomes.r2.1.1.sites.${chr}.vcf.bgz.tbi
    mv gnomad.genomes.r2.1.1.sites.${chr}.vcf.bgz gnomad.genomes.r2.1.1.sites.${chr}.vcf.gz
    mv gnomad.genomes.r2.1.1.sites.${chr}.vcf.bgz.tbi gnomad.genomes.r2.1.1.sites.${chr}.vcf.gz.tbi
done

# Parse INFO columns to retain only relevant ones
bcftools view --header-only gnomad.genomes.r2.1.1.sites.1.vcf.gz | \
    grep -E "AC|AF" | \
    grep -v -e "_male" -e "_female" -e "controls" -e "topmed" -e "neuro" -e "vep" -e "raw" -e "=popmax" -e "AC0" | \
    cut -f3 -d"=" | cut -f1 -d"," \
    > gnomad.genomes.r2.1.1.sites.retained.info

COLS=$(Rscript -e "out = paste0('^INFO/', paste(read.delim('gnomad.genomes.r2.1.1.sites.retained.info', header = F)[['V1']], collapse = ',INFO/')); cat(out)")
CHR=({1..22} X)

# Apply this, and retain filtered sites
for chr in ${CHR[@]}
do
    bcftools annotate \
    --remove "$COLS" \
    --include 'FILTER~"PASS" | FILTER~"RF"' \
    --output-type z \
    --output tmp.${chr}.vcf.gz \
    gnomad.genomes.r2.1.1.sites.${chr}.vcf.gz

    tabix --preset vcf tmp.${chr}.vcf.gz

    bcftools annotate \
    --annotations tmp.${chr}.vcf.gz \
    --include 'FILTER!="PASS"' \
    --mark-sites "+gnomAD_FILTER" \
    -k \
    --output-type z \
    --output gnomad.genomes.r2.1.1.sites.${chr}.minimal.vcf.gz \
    tmp.${chr}.vcf.gz

    tabix --preset vcf gnomad.genomes.r2.1.1.sites.${chr}.minimal.vcf.gz
done

# Concatenate to one file
bcftools concat \
    --output-type z \
    --output gnomad.genomes.r2.1.1.sites.minimal.vcf.gz \
    gnomad.genomes.r2.1.1.sites.{1..22}.minimal.vcf.gz \
    gnomad.genomes.r2.1.1.sites.X.minimal.vcf.gz

tabix --preset vcf gnomad.genomes.r2.1.1.sites.minimal.vcf.gz

# Clean up
for chr in ${CHR[@]}
do
    rm gnomad.genomes.r2.1.1.sites.${chr}.minimal.vcf.gz*
    rm tmp.${chr}.vcf.gz*
done
```