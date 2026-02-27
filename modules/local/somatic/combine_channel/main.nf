process SOMATIC_COMBINE_CHANNEL {
    tag "${meta.id}"
    label 'process_medium'

    conda "bioconda::bcftools=1.9 bioconda::htslib=1.9 bioconda::vt=0.57721 pip::pysam=0.15.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://docker.io/cmopipeline/bcftools-vt:1.2.3' :
        'docker.io/cmopipeline/bcftools-vt:1.2.3' }"

    input:
    tuple val(meta), path(mutect_vcf), path(mutect_tbi),
          path(tumor_bam), path(tumor_bai),
          path(normal_bam), path(normal_bai),
          path(strelka_vcf), path(strelka_tbi)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)
    path repeat_masker
    path repeat_masker_tbi
    path mapability_blacklist
    path mapability_blacklist_tbi
    path pon
    path pon_tbi
    path gnomad
    path gnomad_tbi

    output:
    tuple val(meta), path("${prefix}.pass.vcf"),                            emit: pass_vcf
    tuple val(meta), path("${prefix}.union.annot.vcf.gz"),
                     path("${prefix}.union.annot.vcf.gz.tbi"),              emit: union_annot
    tuple val(meta), path("${prefix}.union.annot.filter.vcf.gz"),
                     path("${prefix}.union.annot.filter.vcf.gz.tbi"),       emit: union_annot_filter
    tuple val(meta), path("${prefix}.union.annot.filter.pass.vcf.gz"),
                     path("${prefix}.union.annot.filter.pass.vcf.gz.tbi"),  emit: union_annot_filter_pass
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def isec_dir = "${meta.tumor_id}.isec"
    """
    # --- VCF header definitions ---
    echo -e '##INFO=<ID=MuTect2,Number=0,Type=Flag,Description="Variant was called by MuTect2">' > vcf.header
    echo -e '##INFO=<ID=Strelka2,Number=0,Type=Flag,Description="Variant was called by Strelka2">' >> vcf.header
    echo -e '##INFO=<ID=Strelka2FILTER,Number=0,Type=Flag,Description="Variant failed filters in Strelka2">' >> vcf.header
    echo -e '##INFO=<ID=RepeatMasker,Number=1,Type=String,Description="RepeatMasker">' > vcf.rm.header
    echo -e '##INFO=<ID=EncodeDacMapability,Number=1,Type=String,Description="EncodeDacMapability">' > vcf.map.header
    echo -e '##INFO=<ID=PoN,Number=1,Type=Integer,Description="Count in panel of normals">' > vcf.pon.header
    echo -e '##FORMAT=<ID=alt_count_raw,Number=1,Type=Integer,Description="Raw alternate allele depth">' > vcf.ad_n.header
    echo -e '##FORMAT=<ID=alt_count_raw_fwd,Number=1,Type=Integer,Description="Raw alternate allele depth on forward strand">' >> vcf.ad_n.header
    echo -e '##FORMAT=<ID=alt_count_raw_rev,Number=1,Type=Integer,Description="Raw alternate allele depth on reverse strand">' >> vcf.ad_n.header
    echo -e '##FORMAT=<ID=ref_count_raw,Number=1,Type=Integer,Description="Raw reference allele depth">' >> vcf.ad_n.header
    echo -e '##FORMAT=<ID=ref_count_raw_fwd,Number=1,Type=Integer,Description="Raw reference allele depth on forward strand">' >> vcf.ad_n.header
    echo -e '##FORMAT=<ID=ref_count_raw_rev,Number=1,Type=Integer,Description="Raw reference allele depth on reverse strand">' >> vcf.ad_n.header
    echo -e '##FORMAT=<ID=depth_raw,Number=1,Type=Integer,Description="Raw total allele depth">' >> vcf.ad_n.header
    echo -e '##FORMAT=<ID=depth_raw_fwd,Number=1,Type=Integer,Description="Raw total allele depth on forward strand">' >> vcf.ad_n.header
    echo -e '##FORMAT=<ID=depth_raw_rev,Number=1,Type=Integer,Description="Raw total allele depth on reverse strand">' >> vcf.ad_n.header

    # --- bcftools isec: get set differences ---
    # 0000: MuTect2 only
    # 0001: Strelka2 only
    # 0002: MuTect2 calls shared by Strelka2
    # 0003: Strelka2 calls shared by MuTect2
    bcftools isec \\
        --output-type z \\
        --prefix ${isec_dir} \\
        ${mutect_vcf} ${strelka_vcf}

    # Annotate shared Strelka2 calls with FILTER status
    bcftools annotate \\
        --annotations ${isec_dir}/0003.vcf.gz \\
        --include 'FILTER!="PASS"' \\
        --mark-sites "+Strelka2FILTER" \\
        -k \\
        --output-type z \\
        --output ${isec_dir}/0003.annot.vcf.gz \\
        ${isec_dir}/0003.vcf.gz

    # Mark MuTect2-only calls
    bcftools annotate \\
        --header-lines vcf.header \\
        --annotations ${isec_dir}/0000.vcf.gz \\
        --mark-sites +MuTect2 \\
        --output-type z \\
        --output ${isec_dir}/0000.annot.vcf.gz \\
        ${isec_dir}/0000.vcf.gz

    # Mark shared calls (MuTect2+Strelka2)
    bcftools annotate \\
        --header-lines vcf.header \\
        --annotations ${isec_dir}/0002.vcf.gz \\
        --mark-sites "+MuTect2;Strelka2" \\
        --output-type z \\
        --output ${isec_dir}/0002.tmp.vcf.gz \\
        ${isec_dir}/0002.vcf.gz

    tabix --preset vcf ${isec_dir}/0002.tmp.vcf.gz
    tabix --preset vcf ${isec_dir}/0003.annot.vcf.gz

    bcftools annotate \\
        --annotations ${isec_dir}/0003.annot.vcf.gz \\
        --columns +INFO,+FORMAT,Strelka2FILTER \\
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
        --output ${meta.tumor_id}.union.vcf.gz

    tabix --preset vcf ${meta.tumor_id}.union.vcf.gz

    # --- gnomAD annotation ---
    bcftools annotate \\
        --annotations ${gnomad} \\
        --columns INFO \\
        --output-type z \\
        --output ${meta.tumor_id}.union.gnomad.vcf.gz \\
        ${meta.tumor_id}.union.vcf.gz

    tabix --preset vcf ${meta.tumor_id}.union.gnomad.vcf.gz

    # --- PoN annotation + flanking sequence ---
    bcftools annotate \\
        --header-lines vcf.pon.header \\
        --annotations ${pon} \\
        --columns PoN:=AC_Het \\
        ${meta.tumor_id}.union.gnomad.vcf.gz | \\
    vt annotate_indels \\
        -r ${fasta} \\
        -o ${meta.tumor_id}.union.annot.vcf -

    mv ${meta.tumor_id}.union.annot.vcf ${prefix}.union.annot.vcf

    # --- Custom filter annotation ---
    filter-vcf.py ${prefix}.union.annot.vcf

    # --- Extract PASS variants ---
    bcftools filter \\
        --include 'FILTER="PASS"' \\
        --output-type v \\
        --output ${prefix}.vcf \\
        ${prefix}.union.annot.filter.vcf

    # --- Re-genotype with GetBaseCountsMultiSample ---
    GetBaseCountsMultiSample \\
        --thread ${task.cpus} \\
        --maq 0 \\
        --fasta ${fasta} \\
        --bam ${meta.tumor_id}:${tumor_bam} \\
        --bam ${meta.normal_id}:${normal_bam} \\
        --vcf ${prefix}.vcf \\
        --output ${prefix}.genotyped.vcf

    bgzip ${prefix}.vcf
    bgzip ${prefix}.genotyped.vcf
    tabix --preset vcf ${prefix}.vcf.gz
    tabix --preset vcf ${prefix}.genotyped.vcf.gz

    # --- Add raw read counts from genotyping ---
    bcftools annotate \\
        --annotations ${prefix}.genotyped.vcf.gz \\
        --header-lines vcf.ad_n.header \\
        --columns FORMAT/alt_count_raw:=FORMAT/AD,FORMAT/ref_count_raw:=FORMAT/RD,FORMAT/alt_count_raw_fwd:=FORMAT/ADP,FORMAT/ref_count_raw_fwd:=FORMAT/RDP,FORMAT/alt_count_raw_rev:=FORMAT/ADN,FORMAT/ref_count_raw_rev:=FORMAT/RDN,FORMAT/depth_raw:=FORMAT/DP,FORMAT/depth_raw_fwd:=FORMAT/DPP,FORMAT/depth_raw_rev:=FORMAT/DPN \\
        --output-type v \\
        --output ${prefix}.union.annot.filter.pass.vcf \\
        ${prefix}.vcf.gz

    cp ${prefix}.union.annot.filter.pass.vcf ${prefix}.pass.vcf

    # --- Compress and index outputs ---
    bgzip ${prefix}.union.annot.vcf
    bgzip ${prefix}.union.annot.filter.vcf
    bgzip ${prefix}.union.annot.filter.pass.vcf
    tabix --preset vcf ${prefix}.union.annot.vcf.gz
    tabix --preset vcf ${prefix}.union.annot.filter.vcf.gz
    tabix --preset vcf ${prefix}.union.annot.filter.pass.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version | head -1 | sed 's/bcftools //')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.pass.vcf
    echo "" | gzip > ${prefix}.union.annot.vcf.gz
    touch ${prefix}.union.annot.vcf.gz.tbi
    echo "" | gzip > ${prefix}.union.annot.filter.vcf.gz
    touch ${prefix}.union.annot.filter.vcf.gz.tbi
    echo "" | gzip > ${prefix}.union.annot.filter.pass.vcf.gz
    touch ${prefix}.union.annot.filter.pass.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: 1.20
    END_VERSIONS
    """
}
