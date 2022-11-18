process GermlineRunSvABA {
    tag "${idNormal}"
    publishDir "${params.outDir}/germline/${idNormal}/svaba", mode: params.publishDirMode, pattern: "*.{vcf.gz,vcf.gz.tbi}"

    input:
    tuple val(idNormal), val(target), path(bamNormal), path(baiNormal), path(targetsBed)
    path(genomeFile)
    path(genomeIndex)
    path(genomeDict)
    path(bwaIndex)

    output:
    tuple val(idNormal), val(target), path("${outputPrefix}.reheader.svaba.germline.sv.vcf.gz"), path("${outputPrefix}.reheader.svaba.germline.sv.vcf.gz.tbi"), emit: SvABA4Combine
    path("*.vcf.gz*"), emit: allVcfs
    path("*.log"), emit: logs
    path("*.txt.gz"), emit: supportingFiles

    script:
    outputPrefix = "${idNormal}"
    target_param = params.assayType == "genome" ? "" : "-k ${targetsBed} "
    """
    svaba run \\
      -t "${bamNormal}" \\
      -G "${genomeFile}" \\
      -p "${task.cpus * 2}" \\
      -I \\
      -L 6 \\
      --id-string "${outputPrefix}" \\
      ${target_param} \\
      -z

    echo -e "${bamNormal} ${idNormal}" > svaba.samplenames.tsv
    bcftools reheader \\
      --samples svaba.samplenames.tsv \\
      --output ${outputPrefix}.reheader.svaba.germline.sv.vcf.gz \\
      ${outputPrefix}.svaba.sv.vcf.gz
    bcftools index -f -t ${outputPrefix}.reheader.svaba.germline.sv.vcf.gz
    """
}
