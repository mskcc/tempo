process SomaticRunSvABA {
	tag "${idTumor}__${idNormal}"
    publishDir "${params.outDir}/somatic/${idTumor}__${idNormal}/svaba", mode: params.publishDirMode, pattern: "*.{vcf.gz,vcf.gz.tbi,contigs.bam}"

    input:
    tuple val(idTumor), val(idNormal), val(target), path(bamTumor), path(baiTumor), path(bamNormal), path(baiNormal)
    path(genomeFile)
    path(genomeIndex)
    path(genomeDict)
    path(bwaIndex)

    output:
    tuple val(idTumor), val(idNormal), val(target), path("*.svaba.somatic.sv.vcf.gz"), path("*.svaba.somatic.sv.vcf.gz.tbi"), emit: SvABA4Combine
    path("*.vcf.gz*"), emit: allVcfs
    path("*.log"), emit: logs
    path("*.txt.gz"), emit: supportingFiles
    path("*.contigs.bam"), emit: contigBams

    script:
    """
    svaba run \\
    -t "${bamTumor}" \\
    -n "${bamNormal}" \\
    -G "${genomeFile}" \\
    -p "${task.cpus * 2}" \\
    --id-string "${idTumor}__${idNormal}" \\
    -z

    rm -f *germline*
    """
}
