process AGGREGATE_SOMATIC_MAF {
    label 'process_single'
    container 'ubuntu:22.04'

    input:
    path(maf_files)

    output:
    path("mut_somatic.maf"), emit: aggregated_maf

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # Extract header from first non-empty MAF file
    head -n 1 \$(for f in ${maf_files}; do [ -s "\$f" ] && echo "\$f" && break; done) > mut_somatic.maf
    
    # Merge all MAF files excluding headers and sort by chromosome and position
    for f in ${maf_files}; do
        tail -n +2 "\$f" 2>/dev/null || true
    done | sort -t'\t' -k5,5 -k6,6n >> mut_somatic.maf
    """

    stub:
    """
    touch mut_somatic.maf
    """
}
