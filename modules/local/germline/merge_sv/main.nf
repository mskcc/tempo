process GERMLINE_MERGE_SV {
    tag "${meta.id}"
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://docker.io/cmopipeline/bcftools-vt-mergesvvcf:0.0.1' :
        'docker.io/cmopipeline/bcftools-vt-mergesvvcf:0.0.1' }"

    input:
    tuple val(meta), path(vcfs), path(tbis), val(caller_names)

    output:
    tuple val(meta), path("${prefix}.merged.vcf.gz"), path("${prefix}.merged.vcf.gz.tbi"), emit: vcf
    path("${prefix}.merged.raw.vcf.{gz,gz.tbi}"), emit: raw_vcf
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"

    // Build label param and input VCFs sorted by caller name (matching original)
    def callers = caller_names instanceof List ? caller_names : [caller_names]
    def vcf_list = vcfs instanceof List ? vcfs : [vcfs]
    def vcf_map = [:]
    for (int i = 0; i < callers.size(); i++) {
        vcf_map[callers[i]] = vcf_list[i]
    }
    def sorted_callers = callers.sort(false)
    def labelparam = sorted_callers.join(",")
    def inVCFs = sorted_callers.collect { vcf_map[it] }.join(" ")
    def passMin = callers.size() > 2 ? 2 : 1
    def sample_id = meta.sample ?: meta.id

    """
    mergesvvcf \\
        -n -m 1 \\
        -l ${labelparam} \\
        -o ${prefix}.merged.raw.vcf \\
        -f -d -s -v \\
        ${inVCFs}

    cat ${prefix}.merged.raw.vcf | \\
        awk -F"\\t" -v OFS="\\t" '\$1 ~ /^#/ && \$1 !~ /^##/ && \$1 !~ /^#CHROM/{next;}{for(i=1; i<=NF; i++) if(\$i ~ /^ *\$/) \$i = "."; print \$0}' | \\
    bcftools sort --temp-dir ./ \\
        > ${prefix}.merged.clean.anon.vcf

    bcftools annotate \\
        --set-id 'TEMPO_%INFO/SVTYPE\\_%CHROM\\_%POS' \\
        -o ${prefix}.merged.clean.vcf \\
        ${prefix}.merged.clean.anon.vcf

    filter-sv-vcf.py \\
        --input ${prefix}.merged.clean.vcf \\
        --output ${prefix}.merged.clean.corrected.vcf \\
        --min ${passMin}

    bcftools view \\
        --samples ${sample_id} \\
        --output-type z \\
        --output-file ${prefix}.merged.vcf.gz \\
        ${prefix}.merged.clean.corrected.vcf

    tabix --preset vcf ${prefix}.merged.vcf.gz

    bcftools view -O z -o ${prefix}.merged.raw.vcf.gz ${prefix}.merged.raw.vcf
    tabix --preset vcf ${prefix}.merged.raw.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version | head -1 | sed 's/bcftools //')
        mergesvvcf: \$(pip show mergesvvcf 2>/dev/null | grep Version | sed 's/Version: //' || echo '1.0.2')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.merged.vcf.gz
    touch ${prefix}.merged.vcf.gz.tbi
    echo "" | gzip > ${prefix}.merged.raw.vcf.gz
    touch ${prefix}.merged.raw.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: 1.20
        mergesvvcf: 1.0.2
    END_VERSIONS
    """
}
