process FACETS_PREVIEW_QC {
    tag "${meta.tumor_id}__${meta.normal_id}"
    label 'process_medium'
    container 'cmopipeline/facets-suite-preview-htstools:0.0.1'

    input:
    tuple val(meta), path(facets_output_files)

    output:
    tuple val(meta), path("*.facets_preview_qc.txt"), emit: facets_qc

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = "${meta.tumor_id}__${meta.normal_id}"
    """
    mkdir -p facets_qc_output

    # Run facetsPreview QC generation
    Rscript - ${args} << 'RSCRIPT'
    library(facetsPreview)
    library(tidyverse)

    output_file <- "${prefix}.facets_preview_qc.txt"

    # Generate genomic annotations and QC metrics
    qc_results <- facetsPreview::generate_genomic_annotations(
        facets_directory = "."
    )

    # Write QC results to file
    write.table(qc_results, output_file, sep = "\\t", row.names = FALSE, quote = FALSE)

    RSCRIPT
    """

    stub:
    def prefix = "${meta.tumor_id}__${meta.normal_id}"
    """
    touch ${prefix}.facets_preview_qc.txt
    """
}
