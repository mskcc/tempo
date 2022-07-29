process SomaticRunLumpy {
	tag "${idTumor}__${idNormal}"
    publishDir "${params.outDir}/somatic/${idTumor}__${idNormal}/lumpy", mode: params.publishDirMode, pattern: "*.lumpy{.vcf.gz,.vcf.gz.tbi}"

    input:
    tuple val(idTumor), val(idNormal), val(target), 
      path(bamTumor), path(baiTumor), 
      path(bamNormal), path(baiNormal)
    path(excludeRegionsBed)
    path(genomeIndex)

    output:
    tuple val(idTumor), val(idNormal), val(target), 
      path("${idTumor}__${idNormal}.lumpy.filtered.vcf.gz"), 
      path("${idTumor}__${idNormal}.lumpy.filtered.vcf.gz.tbi")
    
    script:
    """
    echo -e "${idTumor}\t${bamTumor}\n${idNormal}\t${bamNormal}" > inputs.txt
    cat inputs.txt | while read label bam ; do 
    samtools view \\
        -@ ${task.cpus} \\
        -b -F 1294 \\
        \${bam} > \\
        \${label}.discordants.unsorted.bam
    samtools sort \\
        -@ ${task.cpus} \\
        -o \${label}.discordants.bam \\
        \${label}.discordants.unsorted.bam 

    samtools view \\
        -@ ${task.cpus} \\
        -h \${bam} \\
        | /opt/lumpy-sv/scripts/extractSplitReads_BwaMem -i stdin \\
        | samtools view -Sb - \\
        > \${label}.splitters.unsorted.bam
    samtools sort \\
        -@ ${task.cpus} \\
        -o \${label}.splitters.bam \\
        \${label}.splitters.unsorted.bam 
    done
    rm -f inputs.txt

    lumpyexpress \\
	  -B ${bamTumor},${bamNormal} \\
	  -S ${idTumor}.splitters.bam,${idNormal}.splitters.bam \\
      -D ${idTumor}.discordants.bam,${idNormal}.discordants.bam \\
      -x ${excludeRegionsBed} \\
	  -o ${idTumor}__${idNormal}.raw.vcf

    svtyper \\
      -i ${idTumor}__${idNormal}.raw.vcf \\
      -B ${bamTumor},${bamNormal} \\
      -l stats.json > \\
      ${idTumor}__${idNormal}.gt.vcf

    awk '{printf "##contig=<ID=%s,length=%s>\\n", \$1, \$2}' ${genomeIndex} > contigs.txt

    cat \\
      <(head -3 ${idTumor}__${idNormal}.gt.vcf) \\
      contigs.txt \\
      <(tail -n +4 ${idTumor}__${idNormal}.gt.vcf) > \\
      ${idTumor}__${idNormal}.lumpy.vcf

    bcftools filter \\
      -i "FORMAT/AB[1:0] < .25 | FORMAT/AB[1:0]=\\".\\"" \\
      ${idTumor}__${idNormal}.lumpy.vcf | \\
    bcftools filter \\
      -i "FORMAT/AB[0:0] >= .05" | \\
    bcftools filter \\
      -i "QUAL>20" > \\
      ${idTumor}__${idNormal}.lumpy.filtered.vcf

    bcftools sort ${idTumor}__${idNormal}.lumpy.filtered.vcf -Oz -o ${idTumor}__${idNormal}.lumpy.filtered.vcf.gz
    bcftools index -t ${idTumor}__${idNormal}.lumpy.filtered.vcf.gz

    rm -f ${idTumor}__${idNormal}.lumpy.vcf ${idTumor}__${idNormal}.gt.vcf ${idTumor}__${idNormal}.raw.vcf

    """
}
