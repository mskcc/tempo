process GERMLINE_COMBINE_CHANNEL {
    tag "${meta.id}"
    label 'process_medium'

    conda "bioconda::bcftools=1.9 bioconda::htslib=1.9 bioconda::vt=0.57721 pip::pysam=0.15.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://cmopipeline/bcftools-vt:1.2.2' :
        'cmopipeline/bcftools-vt:1.2.2' }"

    input:
    tuple val(meta), path(hc_vcf), path(hc_tbi),
          path(strelka_vcf), path(strelka_tbi),
          path(tumor_bam), path(tumor_bai)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)
    path repeat_masker
    path repeat_masker_tbi
    path mapability_blacklist
    path mapability_blacklist_tbi
    path gnomad
    path gnomad_tbi

    output:
    tuple val(meta), path("${prefix}.germline.vcf"),   emit: germline_vcf
    tuple val(meta), path("${meta.normal_id}.union.vcf.gz"),
                     path("${meta.normal_id}.union.vcf.gz.tbi"),      emit: union
    tuple val(meta), path("${meta.normal_id}.union.pass.vcf.gz"),
                     path("${meta.normal_id}.union.pass.vcf.gz.tbi"), emit: union_pass
    tuple val(meta), path("${prefix}.germline.vcf.gz"),
                     path("${prefix}.germline.vcf.gz.tbi"),           emit: germline_vcf_gz
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def isec_dir = "${meta.normal_id}.isec"
    def gnomad_af_filter = params.germline_gnomad_af ?: 'INFO/AF > 0.0001'
    """
    # --- VCF header definitions ---
    echo -e '##INFO=<ID=HaplotypeCaller,Number=0,Type=Flag,Description="Variant was called by HaplotypeCaller">' > vcf.header
    echo -e '##INFO=<ID=Strelka2,Number=0,Type=Flag,Description="Variant was called by Strelka2">' >> vcf.header
    echo -e '##INFO=<ID=Strelka2FILTER,Number=0,Type=Flag,Description="Variant failed filters in Strelka2">' >> vcf.header
    echo -e '##INFO=<ID=RepeatMasker,Number=1,Type=String,Description="RepeatMasker">' > vcf.rm.header
    echo -e '##INFO=<ID=EncodeDacMapability,Number=1,Type=String,Description="EncodeDacMapability">' > vcf.map.header

    # --- bcftools isec: get set differences ---
    bcftools isec \\
        --output-type z \\
        --prefix ${isec_dir} \\
        ${hc_vcf} ${strelka_vcf}

    # Annotate shared Strelka2 calls with FILTER status
    bcftools annotate \\
        --annotations ${isec_dir}/0003.vcf.gz \\
        --include 'FILTER!="PASS"' \\
        --mark-sites "+Strelka2FILTER" \\
        -k \\
        --output-type z \\
        --output ${isec_dir}/0003.annot.vcf.gz \\
        ${isec_dir}/0003.vcf.gz

    # Mark HaplotypeCaller-only calls
    bcftools annotate \\
        --header-lines vcf.header \\
        --annotations ${isec_dir}/0000.vcf.gz \\
        --mark-sites +HaplotypeCaller \\
        --output-type z \\
        --output ${isec_dir}/0000.annot.vcf.gz \\
        ${isec_dir}/0000.vcf.gz

    # Mark shared calls (HaplotypeCaller+Strelka2)
    bcftools annotate \\
        --header-lines vcf.header \\
        --annotations ${isec_dir}/0002.vcf.gz \\
        --mark-sites "+HaplotypeCaller;Strelka2" \\
        --output-type z \\
        --output ${isec_dir}/0002.tmp.vcf.gz \\
        ${isec_dir}/0002.vcf.gz

    tabix --preset vcf ${isec_dir}/0002.tmp.vcf.gz
    tabix --preset vcf ${isec_dir}/0003.annot.vcf.gz

    bcftools annotate \\
        --annotations ${isec_dir}/0003.annot.vcf.gz \\
        --columns +FORMAT,Strelka2FILTER \\
        --output-type z \\
        --output ${isec_dir}/0002.annot.vcf.gz \\
        ${isec_dir}/0002.tmp.vcf.gz

    # Mark Strelka2-only calls
    bcftools annotate \\
        --header-lines vcf.header \\
        --annotations ${isec_dir}/0001.vcf.gz \\
        --mark-sites +Strelka2 \\
        --output-type z \\
        --output ${isec_dir}/0001.annot.vcf.gz \\
        ${isec_dir}/0001.vcf.gz

    tabix --preset vcf ${isec_dir}/0000.annot.vcf.gz
    tabix --preset vcf ${isec_dir}/0001.annot.vcf.gz
    tabix --preset vcf ${isec_dir}/0002.annot.vcf.gz

    # --- Concatenate and annotate with blacklists ---
    bcftools concat \\
        --allow-overlaps \\
        --rm-dups all \\
        ${isec_dir}/0000.annot.vcf.gz \\
        ${isec_dir}/0001.annot.vcf.gz \\
        ${isec_dir}/0002.annot.vcf.gz | \\
    bcftools sort | \\
    bcftools annotate \\
        --header-lines vcf.rm.header \\
        --annotations ${repeat_masker} \\
        --columns CHROM,FROM,TO,RepeatMasker | \\
    bcftools annotate \\
        --header-lines vcf.map.header \\
        --annotations ${mapability_blacklist} \\
        --columns CHROM,FROM,TO,EncodeDacMapability \\
        --output-type z \\
        --output ${meta.normal_id}.union.vcf.gz

    tabix --preset vcf ${meta.normal_id}.union.vcf.gz

    # --- Filter PASS variants ---
    bcftools filter \\
        --include 'FILTER="PASS"' \\
        --output-type z \\
        --output ${meta.normal_id}.union.pass.vcf.gz \\
        ${meta.normal_id}.union.vcf.gz

    tabix --preset vcf ${meta.normal_id}.union.pass.vcf.gz

    # --- gnomAD annotation + AF filter ---
    bcftools annotate \\
        --annotations ${gnomad} \\
        --columns INFO \\
        ${meta.normal_id}.union.pass.vcf.gz | \\
    bcftools filter \\
        --exclude "${gnomad_af_filter}" \\
        --output-type v \\
        --output ${meta.normal_id}.union.gnomad.vcf

    # --- Re-genotype with tumor BAM ---
    GetBaseCountsMultiSample \\
        --fasta ${fasta} \\
        --bam ${meta.tumor_id}:${tumor_bam} \\
        --vcf ${meta.normal_id}.union.gnomad.vcf \\
        --output ${meta.tumor_id}.genotyped.vcf

    bgzip ${meta.normal_id}.union.gnomad.vcf
    bgzip ${meta.tumor_id}.genotyped.vcf
    tabix --preset vcf ${meta.normal_id}.union.gnomad.vcf.gz
    tabix --preset vcf ${meta.tumor_id}.genotyped.vcf.gz

    # --- Merge normal variants with tumor genotypes ---
    bcftools merge \\
        --output ${prefix}.germline.vcf \\
        --output-type v \\
        ${meta.normal_id}.union.gnomad.vcf.gz \\
        ${meta.tumor_id}.genotyped.vcf.gz

    bgzip -c ${prefix}.germline.vcf > ${prefix}.germline.vcf.gz
    tabix --preset vcf ${prefix}.germline.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version | head -1 | sed 's/bcftools //')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.germline.vcf
    echo "" | gzip > ${meta.normal_id}.union.vcf.gz
    touch ${meta.normal_id}.union.vcf.gz.tbi
    echo "" | gzip > ${meta.normal_id}.union.pass.vcf.gz
    touch ${meta.normal_id}.union.pass.vcf.gz.tbi
    echo "" | gzip > ${prefix}.germline.vcf.gz
    touch ${prefix}.germline.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: 1.20
    END_VERSIONS
    """
}
