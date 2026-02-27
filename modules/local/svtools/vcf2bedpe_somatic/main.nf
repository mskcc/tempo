process SVTOOLS_VCF2BEDPE_SOMATIC {
    tag "${meta.id}"
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://docker.io/cmopipeline/svtools:0.0.3' :
        'docker.io/cmopipeline/svtools:0.0.3' }"

    input:
    tuple val(meta), path(vcf), path(tbi)

    output:
    tuple val(meta), path("${prefix}.combined.bedpe"), emit: bedpe
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    export LC_ALL=C

    echo -e "${meta.tumor_id} TUMOR\\n${meta.normal_id} NORMAL" > normalize.samplenames.tsv
    bcftools reheader \\
        --samples normalize.samplenames.tsv \\
        --output reheader_input.vcf.gz \\
        ${vcf}

    svtools vcftobedpe \\
        -i reheader_input.vcf.gz \\
        -o ${prefix}.combined.tmp.bedpe \\
        -t ${prefix}_tmp

    if [ ! -s ${prefix}.combined.tmp.bedpe ] ; then
        echo -e "#CHROM_A\\tSTART_A\\tEND_A\\tCHROM_B\\tSTART_B\\tEND_B\\tID\\tQUAL\\tSTRAND_A\\tSTRAND_B\\tTYPE\\tFILTER\\tNAME_A\\tREF_A\\tALT_A\\tNAME_B\\tREF_B\\tALT_B\\tINFO_A\\tINFO_B\\tFORMAT\\tTUMOR\\tNORMAL" >> ${prefix}.combined.tmp.bedpe
    fi

    zgrep "^##" ${vcf} | sed "s/##fileformat=*/##fileformat=BEDPE/g" > ${prefix}.combined.unsorted.bedpe
    grep -v "^##" ${prefix}.combined.tmp.bedpe | \\
        awk -F"\\t" -v tid="${meta.tumor_id}" -v nid="${meta.normal_id}" -v OFS="\\t" \\
        'NR == 1 {print \$0,"TUMOR_ID","NORMAL_ID";next;}{print \$0,tid,nid}' \\
        >> ${prefix}.combined.unsorted.bedpe

    svtools bedpesort \\
        ${prefix}.combined.unsorted.bedpe \\
        ${prefix}.combined.bedpe

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        svtools: \$(svtools --version 2>&1 | head -1 || echo "0.0.3")
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.combined.bedpe
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        svtools: 0.0.3
    END_VERSIONS
    """
}
