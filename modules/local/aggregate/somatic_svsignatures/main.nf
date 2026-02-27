// Aggregate Somatic SV Signatures — matches original Tempo SomaticAggregateSvSignatures
// Merges exposure TSVs and concatenates catalogue PDFs via ghostscript
process AGGREGATE_SOMATIC_SVSIGNATURES {
    tag "${cohort}"
    label 'process_single'
    container 'docker.io/cmopipeline/signaturetoolslib:0.0.1'

    input:
    val(cohort)
    path(catalogue_pdfs)
    path(exposure_files)

    output:
    path("sv_catalogues.pdf"), emit: catalogues
    path("sv_exposures.tsv"), emit: exposures

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    awk 'FNR==1 && NR!=1{next;}{print}' *_exposures.tsv > sv_exposures.tsv
    gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -dPDFSETTINGS=/prepress -sOutputFile=new.pdf *_catalogues.pdf
    mv new.pdf sv_catalogues.pdf
    """

    stub:
    """
    touch sv_catalogues.pdf sv_exposures.tsv
    """
}
