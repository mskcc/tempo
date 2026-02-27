// Polysolver HLA typing — matches original Tempo RunPolysolver
// Genome-aware mapping: GRCh37→hg19, GRCh38→hg38, smallGRCh37→hardcoded
process POLYSOLVER {
    tag "$meta.id"
    label 'process_high'

    container "docker.io/sachet/polysolver:v4"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("*.hla.txt"),    emit: hla_types
    tuple val(meta), path("*"),            emit: results
    path "versions.yml",                   emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def genome = params.genome ?: 'GRCh37'
    def genome_mapped = genome == 'GRCh37' ? 'hg19' : (genome == 'GRCh38' ? 'hg38' : (genome == 'smallGRCh37' ? 'small' : genome))
    """
    if [ "${genome_mapped}" != "small" ] ; then

        cp /home/polysolver/scripts/shell_call_hla_type .
        sed -i "171s|TMP_DIR=.*|TMP_DIR=\$(pwd)/nf-scratch/|" shell_call_hla_type

        mkdir -p nf-scratch

        bash shell_call_hla_type \\
            ${bam} \\
            Unknown \\
            1 \\
            ${genome_mapped} \\
            STDFQ \\
            0 \\
            .

        mv winners.hla.txt ${prefix}.hla.txt

    else

        echo -e 'HLA-A\\thla_a_01_01_01_01\\thla_a_01_01_01_01\\nHLA-B\\thla_b_15_02_01\\thla_b_15_02_01\\nHLA-C\\thla_c_01_02_01\\thla_c_01_02_01' > ${prefix}.hla.txt

    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        polysolver: "v4"
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo -e 'HLA-A\\thla_a_01_01_01_01\\thla_a_01_01_01_01\\nHLA-B\\thla_b_15_02_01\\thla_b_15_02_01\\nHLA-C\\thla_c_01_02_01\\thla_c_01_02_01' > ${prefix}.hla.txt
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        polysolver: "stub"
    END_VERSIONS
    """
}
