// Extract read group ID from FASTQ header
// Matches original Tempo AlignReads behavior:
//   rgID=`zcat $fastqFile1 | head -1 | tr ':/' '@' | cut -d '@' -f2-5`
//   readGroup="@RG\\tID:${rgID}\\tSM:${idSample}\\tLB:${idSample}\\tPL:Illumina"
process EXTRACT_READ_GROUP {
    tag "${meta.id}"
    label 'process_single'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), stdout, emit: read_group

    when:
    task.ext.when == null || task.ext.when

    script:
    def fastq = reads instanceof List ? reads[0] : reads
    """
    printf "\$(zcat ${fastq} | head -1 | tr ':/' '@' | cut -d '@' -f2-5)"
    """

    stub:
    """
    printf "${meta.sample}_${meta.lane ?: 'L001'}"
    """
}
