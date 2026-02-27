// LOHHLA — matches original Tempo RunLOHHLA
// Includes HLA file preprocessing, purity/ploidy extraction, output renaming
process LOHHLA {
    tag "$meta.id"
    label 'process_high'

    container "docker.io/cmopipeline/lohhla:1.1.7"

    input:
    tuple val(meta), path(tumor_bam), path(tumor_bai), path(normal_bam), path(normal_bai), path(hla_types), path(purity_out)
    path hla_fasta
    path hla_dat

    output:
    tuple val(meta), path("*.DNA.HLAlossPrediction_CI.txt"), emit: predictions, optional: true
    tuple val(meta), path("*.DNA.IntegerCPN_CI.txt"), emit: integer_cpn, optional: true
    tuple val(meta), path("*.HLA.pdf"), emit: figures, optional: true
    tuple val(meta), path("*"), emit: results
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.tumor_id}__${meta.normal_id}"
    def min_coverage = params.lohhla_min_coverage_filter ?: '10'
    """
    # Massage HLA file: convert tabs to newlines and remove HLA- prefix lines
    cat ${hla_types} | tr "\\t" "\\n" | grep -v "HLA" > massaged.winners.hla.txt

    # Extract purity and ploidy from FACETS output
    PURITY=\$(grep Purity ${purity_out} | grep -oP "[0-9\\.]+|NA+" || echo "NA")
    PLOIDY=\$(grep Ploidy ${purity_out} | grep -oP "[0-9\\.]+|NA+" || echo "NA")
    echo -e "tumorPurity\\ttumorPloidy" > tumor_purity_ploidy.txt
    echo -e "${meta.tumor_id}\\t\$PURITY\\t\$PLOIDY" >> tumor_purity_ploidy.txt

    Rscript --no-init-file /lohhla/LOHHLAscript.R \\
        --patientId ${prefix} \\
        --normalBAMfile ${normal_bam} \\
        --tumorBAMfile ${tumor_bam} \\
        --HLAfastaLoc ${hla_fasta} \\
        --HLAexonLoc ${hla_dat} \\
        --CopyNumLoc tumor_purity_ploidy.txt \\
        --minCoverageFilter ${min_coverage} \\
        --hlaPath massaged.winners.hla.txt \\
        --gatkDir /picard-tools \\
        --novoDir /opt/conda/bin \\
        ${args}

    # Post-processing: rename and fix output files (matches original Tempo)
    if [[ -f ${prefix}.${min_coverage}.DNA.HLAlossPrediction_CI.txt ]] ; then
        sed -i "s/^${meta.tumor_id}/${prefix}/g" ${prefix}.${min_coverage}.DNA.HLAlossPrediction_CI.txt
    else
        rm -rf *.DNA.HLAlossPrediction_CI.txt
    fi

    if [[ -f ${prefix}.${min_coverage}.DNA.IntegerCPN_CI.txt ]] ; then
        sed -i "s/^/${prefix}\\t/g" ${prefix}.${min_coverage}.DNA.IntegerCPN_CI.txt
        sed -i "0,/^${prefix}\\t/s//sample\\t/" ${prefix}.${min_coverage}.DNA.IntegerCPN_CI.txt
    else
        rm -rf *.DNA.IntegerCPN_CI.txt
    fi

    touch ${prefix}.${min_coverage}.DNA.HLAlossPrediction_CI.txt
    touch ${prefix}.${min_coverage}.DNA.IntegerCPN_CI.txt

    mv ${prefix}.${min_coverage}.DNA.HLAlossPrediction_CI.txt ${prefix}.DNA.HLAlossPrediction_CI.txt
    mv ${prefix}.${min_coverage}.DNA.IntegerCPN_CI.txt ${prefix}.DNA.IntegerCPN_CI.txt

    if find Figures -mindepth 1 2>/dev/null | read ; then
        mv Figures/* .
        if [[ -f ${meta.tumor_id}.minCoverage_${min_coverage}.HLA.pdf ]] ; then
            mv ${meta.tumor_id}.minCoverage_${min_coverage}.HLA.pdf ${prefix}.HLA.pdf
        fi
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        lohhla: "1.1.7"
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.tumor_id}__${meta.normal_id}"
    """
    touch ${prefix}.DNA.HLAlossPrediction_CI.txt
    touch ${prefix}.DNA.IntegerCPN_CI.txt
    touch ${prefix}.HLA.pdf
    touch ${prefix}_result.txt
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        lohhla: "stub"
    END_VERSIONS
    """
}
