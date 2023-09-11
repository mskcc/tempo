process SomaticRunSvABA {
	tag "${idTumor}__${idNormal}"
    publishDir "${params.outDir}/somatic/${idTumor}__${idNormal}/svaba", mode: params.publishDirMode, pattern: "*.{vcf.gz,vcf.gz.tbi}"

    input:
    tuple val(idTumor), val(idNormal), val(target), path(bamTumor), path(baiTumor), path(bamNormal), path(baiNormal), path(targetsBed)
    path(genomeFile)
    path(genomeIndex)
    path(genomeDict)
    path(bwaIndex)

    output:
    tuple val(idTumor), val(idNormal), val(target), path("${outputPrefix}.reheader.svaba.somatic.sv.vcf.gz"), path("${outputPrefix}.reheader.svaba.somatic.sv.vcf.gz.tbi"), emit: SvABA4Combine
    tuple val(idTumor), val(idNormal), val(target), path("${outputPrefix}.INDEL*vcf.gz"), path("${outputPrefix}.INDEL*vcf.gz.tbi"), emit: SvABA4INDELCombine
    path("*.vcf.gz*"), emit: allVcfs
    path("*.log"), emit: logs
    path("*.txt.gz"), emit: supportingFiles

    script:
    outputPrefix = "${idTumor}__${idNormal}"
    target_param = params.assayType == "genome" ? "" : "-k ${targetsBed} "
    """
    svaba run \\
      -t "${bamTumor}" \\
      -n "${bamNormal}" \\
      -G "${genomeFile}" \\
      -p "${task.cpus * 2}" \\
      --id-string "${outputPrefix}" \\
      ${target_param} \\
      -z

    rm -f *germline*

    echo -e "${bamTumor} ${idTumor}\\n${bamNormal} ${idNormal}" > svaba.samplenames.tsv
    bcftools reheader \\
      --samples svaba.samplenames.tsv \\
      --output ${outputPrefix}.reheader.svaba.somatic.sv.vcf.gz \\
      ${outputPrefix}.svaba.somatic.sv.vcf.gz
    bcftools index -f -t ${outputPrefix}.reheader.svaba.somatic.sv.vcf.gz
    """
}
