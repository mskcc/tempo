nextflow.enable.dsl=2

params.tumor_bam = 'tumor.bam'
params.tumor_bai = 'tumor.bam.bai'
params.normal_bam = 'normal.bam'
params.normal_bai = 'normal.bam.bai'

params.genome_base =  "/juno/work/taylorlab/cmopipeline/mskcc-igenomes/igenomes/Homo_sapiens/GATK/GRCh37"
bwaIndex         = "${params.genome_base}/Sequence/BWAIndex/human_g1k_v37_decoy.fasta.{amb,ann,bwt,pac,sa}"
genomeDict       = "${params.genome_base}/Sequence/WholeGenomeFasta/human_g1k_v37_decoy.dict"
genomeFile       = "${params.genome_base}/Sequence/WholeGenomeFasta/human_g1k_v37_decoy.fasta"
genomeIndex      = "${params.genome_base}/Sequence/WholeGenomeFasta/human_g1k_v37_decoy.fasta.fai"

process run_svaba {
    // https://github.com/walaj/svaba#whole-genome-somatic-sv-and-indel-detection
    input:
    tuple path(tumor_bam), path(tumor_bai), path(normal_bam), path(normal_bai), path(genome_fasta), file(genome_files: "*"), file(bwa_files: "*")

    script:
    dedup_tumor_bam = "${tumor_bam}.dedup.bam"
    dedup_normal_bam = "${normal_bam}.dedup.bam"
    """
    svaba --version

    # samtools 1.2
    # Using htslib 1.2.1

    echo "samtools view on "${tumor_bam}""
    date
    samtools view -F 1024 -b -o "${dedup_tumor_bam}" "${tumor_bam}"

    echo "samtools view on "${normal_bam}""
    date
    samtools view -F 1024 -b -o "${dedup_normal_bam}" "${normal_bam}"

    echo "samtools index on "${dedup_tumor_bam}""
    date
    samtools index "${dedup_tumor_bam}"

    echo "samtools index on "${dedup_normal_bam}""
    date
    samtools index "${dedup_normal_bam}"

    echo "run svaba"
    date
    svaba run \
    -t "${dedup_tumor_bam}" \
    -n "${dedup_normal_bam}" \
    -G "${genome_fasta}" \
    -p "${task.cpus}" \
    -z

    echo "svaba done"
    date
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

workflow {
    tumor_bam = Channel.fromPath("${params.tumor_bam}")
    tumor_bai = Channel.fromPath("${params.tumor_bai}")
    normal_bam = Channel.fromPath("${params.normal_bam}")
    normal_bai = Channel.fromPath("${params.normal_bai}")

    tumor_bam.combine(tumor_bai).combine(normal_bam).combine(normal_bai).set{ sample_files }

    genome_fasta = Channel.fromPath("${genomeFile}")

    genome_files = Channel.from([
        file("${genomeIndex}"),
        file("${genomeDict}"),
        ]).collect().toList()

    bwa_indexes = Channel.fromPath("${bwaIndex}").collect().toList()

    sample_files.combine(genome_fasta).combine(genome_files).combine(bwa_indexes).set { input_files }

    // input_files.view()
    run_svaba(input_files)
}
