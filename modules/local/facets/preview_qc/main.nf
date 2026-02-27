// FACETS Preview QC — matches original Tempo DoFacetsPreviewQC
// Uses facetsPreview::generate_genomic_annotations with config file
process FACETS_PREVIEW_QC {
    tag "${meta.tumor_id}__${meta.normal_id}"
    label 'process_medium'
    container 'docker.io/cmopipeline/facets-suite-preview-htstools:0.0.1'

    input:
    tuple val(meta), path(facets_output_files), path(counts_file)

    output:
    tuple val(meta), path("*.facets_preview_qc.txt"), emit: facets_qc

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = "${meta.tumor_id}__${meta.normal_id}"
    """
    mkdir -p facets_output
    for i in ${facets_output_files} ; do
        cp \$i facets_output/\$i
    done

    echo -e "sample_id\\tsample_path\\ttumor_id" > manifest.txt
    echo -e "${prefix}\\t\$(pwd)\\t${meta.tumor_id}" >> manifest.txt
    gzip manifest.txt

    mkdir -p refit_watcher/bin/ refit_watcher/refit_jobs/

    R -e "facetsPreview::generate_genomic_annotations('${prefix}', '\$(pwd)/', '/usr/bin/facets-preview/tempo_config.json')"

    cp facets_qc.txt ${prefix}.facets_preview_qc.txt
    rm -f facets_output/*
    """

    stub:
    def prefix = "${meta.tumor_id}__${meta.normal_id}"
    """
    touch ${prefix}.facets_preview_qc.txt
    """
}
