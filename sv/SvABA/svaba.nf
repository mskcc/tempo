nextflow.enable.dsl=2


process run_svaba {
    // run SvABA
    // https://github.com/walaj/svaba#whole-genome-somatic-sv-and-indel-detection
    input:
    tuple path(tumor_bam), path(tumor_bai), path(normal_bam), path(normal_bai), path(genome_fasta), file(genome_files: "*"), file(bwa_files: "*")

    // TODO: what to output?

    script:
    """
    svaba --version
    svaba run \
    -t "${tumor_bam}" \
    -n "${normal_bam}" \
    -G "${genome_fasta}" \
    -p "${task.cpus}" \
    -z
    """

    // output:
    // no_id.log
    // no_id.alignments.txt.gz
    // no_id.bps.txt.gz
    // no_id.discordant.txt.gz
    // no_id.svaba.unfiltered.germline.indel.vcf.gz
    // no_id.svaba.unfiltered.somatic.indel.vcf.gz
    // no_id.svaba.unfiltered.germline.sv.vcf.gz
    // no_id.svaba.unfiltered.somatic.sv.vcf.gz
    // no_id.svaba.germline.indel.vcf.gz
    // no_id.svaba.somatic.indel.vcf.gz
    // no_id.svaba.germline.sv.vcf.gz
    // no_id.svaba.somatic.sv.vcf.gz
    // no_id.contigs.bam
}

workflow svaba_workflow {
    // sub-workflow unit to run SvABA on a single tumor normal pair
    take:
        sample_files // [tumor_bam, tumor_bai, normal_bam, normal_bai]
        genome_fasta
        genome_files // [genomeIndex, genomeDict]
        bwa_indexes // [.fasta.{amb,ann,bwt,pac,sa}]

    main:
        sample_files.combine(genome_fasta).combine(genome_files).combine(bwa_indexes).set { input_files }
        run_svaba(input_files)

    // emit:
    // TODO: what to output?
}


workflow {
    // default workflow when this script is run from command line
    params.tumor_bam = 'tumor.bam'
    params.tumor_bai = 'tumor.bam.bai'
    params.normal_bam = 'normal.bam'
    params.normal_bai = 'normal.bam.bai'

    tumor_bam = Channel.fromPath("${params.tumor_bam}")
    tumor_bai = Channel.fromPath("${params.tumor_bai}")
    normal_bam = Channel.fromPath("${params.normal_bam}")
    normal_bai = Channel.fromPath("${params.normal_bai}")

    tumor_bam.combine(tumor_bai).combine(normal_bam).combine(normal_bai).set{ sample_files }
    genome_fasta = Channel.fromPath("${params.genomeFile}")
    genome_files = Channel.from([ file("${params.genomeIndex}"), file("${params.genomeDict}") ]).collect().toList()
    bwa_indexes = Channel.fromPath("${params.bwaIndex}").collect().toList()

    svaba_workflow(sample_files, genome_fasta, genome_files, bwa_indexes)
}
