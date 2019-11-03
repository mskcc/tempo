# Gene Annotation for SVs

Gene annotation from [Gencode](https://www.gencodegenes.org) is used for functional-effect prediction of SVs. This annotation is retrieved and parsed as such:

```shell
# Download
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_31/GRCh37_mapping/gencode.v31lift37.basic.annotation.gtf.gz

# Subset on protein-coding genes, remove "chr" prefix
# Also remove IGKV1D which for some reason causes a bug in svtk
zgrep "protein_coding" gencode.v31lift37.basic.annotation.gtf.gz | \
    sed 's/^chr//g' | \
    grep -v "IGKV1D" \
    > gencode.v31lift37.basic.annotation.protein_coding.gtf
```
