process BWAMEM2_INDEX {
    tag "${meta.id}"
    label 'process_high'

    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e0/e05ce34b46ad42810eb29f74e4e304c0cb592b2ca15572929ed8bbaee58faf01/data'
        : 'community.wave.seqera.io/library/bwa-mem2_htslib_samtools:db98f81f55b64113'}"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("bwamem2"), emit: index
    path "versions.yml",              emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    mkdir bwamem2
    bwa-mem2 index ${fasta} -p bwamem2/${fasta.getName()}
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwamem2: \$(bwa-mem2 version 2>/dev/null || echo "N/A")
    END_VERSIONS
    """

    stub:
    """
    mkdir bwamem2
    touch bwamem2/${fasta.getName()}.amb
    touch bwamem2/${fasta.getName()}.ann
    touch bwamem2/${fasta.getName()}.bwt.2bit.64
    touch bwamem2/${fasta.getName()}.pac
    touch bwamem2/${fasta.getName()}.0123
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwamem2: stub
    END_VERSIONS
    """
}
