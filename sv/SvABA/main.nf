nextflow.enable.dsl=2

include { svaba_workflow } from './svaba.nf'

workflow {
    // workflow to run SvABA on a batch of tumor normal pairs
    genome_fasta = Channel.fromPath("${params.genomeFile}")
    genome_files = Channel.from([ file("${params.genomeIndex}"), file("${params.genomeDict}") ]).collect().toList()
    bwa_indexes = Channel.fromPath("${params.bwaIndex}").collect().toList()

    sample_files = Channel.from([
        [ file('input/tumor.bam'), file('input/tumor.bam.bai'), file('input/normal.bam'), file('input/normal.bam.bai') ],
        [ file('input/tumor2.bam'), file('input/tumor2.bam.bai'), file('input/normal2.bam'), file('input/normal2.bam.bai') ],
        [ file('input/tumor3.bam'), file('input/tumor3.bam.bai'), file('input/normal3.bam'), file('input/normal3.bam.bai') ]
        ])

    svaba_workflow(sample_files, genome_fasta, genome_files, bwa_indexes)
}
