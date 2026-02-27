# Pipeline Output

For a detailed description of the pipeline output files and their formats, please see the [Output Documentation](./outputs.md).

The TEMPO pipeline generates comprehensive genomic analysis results including:

- **Alignment results**: BAM files with aligned reads and quality metrics
- **Variant calls**: VCF files for somatic and germline variants
- **Copy number analysis**: Segmentation files from FACETS
- **Structural variants**: VCF files for detected structural variants
- **Quality control**: MultiQC reports summarizing analysis quality
- **Annotation**: Annotated VCF and MAF files with functional predictions
- **Immunogenomics**: HLA typing and neoantigen predictions

All output files are organized in the specified output directory with consistent naming conventions for easy downstream analysis.
