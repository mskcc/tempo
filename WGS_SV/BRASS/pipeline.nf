genome_base      = params.genome_base
genomeFile       = "${params.genome_base}/Sequence/WholeGenomeFasta/human_g1k_v37_decoy.fasta"
genomeIndex      = "${params.genome_base}/Sequence/WholeGenomeFasta/human_g1k_v37_decoy.fasta.fai"
genomeDict       = "${params.genome_base}/Sequence/WholeGenomeFasta/human_g1k_v37_decoy.dict"
bwaIndex         = "${params.genome_base}/Sequence/BWAIndex/human_g1k_v37_decoy.fasta.{amb,ann,bwt,pac,sa}"
tumor_path       = "tumor.bam"
normal_path      = "normal.bam"
RefBase          = "/juno/work/ccs/noronhaa/BRASS_test/ref"

Channel.value([file(tumor_path), file(tumor_path + ".bai"), file(normal_path), file(normal_path + ".bai")]).view().into{bamInput; bamInput2}
Channel.from([file(tumor_path), file(tumor_path + ".bai")], [file(normal_path), file(normal_path + ".bai")]).set{bamInput3}

process runPCAP {
	
	input:
	set file(bam), file(bai) from bamInput3
	set file(genomeFile), file(genomeIndex), file(bwaIndexFiles), file(genomeDict), file(REF_BASE) from Channel.value([file(genomeFile), file(genomeIndex), file(bwaIndex), file(genomeDict), file(RefBase)])

	output: 
	file("*.bas") into basFiles
	file("listoffiles") into pcapOutFiles

	script:
	"""
	mkdir myout
	bam_stats -i ${bam} -o ${bam}.bas -r ${genomeIndex} -@ ${ task.cpus > 1 ? task.cpus - 1 : task.cpus }
	find . > listoffiles
	"""
}

basFiles.collect().set{basFiles}

process runAscat {
	label 'ascat' 
	echo true

	input: 
	set file(tumorBam), file(tumorBai), file(normalBam), file(normalBai) from bamInput
	set file(genomeFile), file(genomeIndex), file(REF_BASE) from Channel.value([file(genomeFile), file(genomeIndex), file(RefBase)])

	output:
	file("ascatResults/*samplestatistics.txt") into ascatOut4Brass
	file("listoffiles") into ascatOutFiles

	script:
	// -l $REF_BASE/gender.tsv       in its absence refers to share/gender/GRCh37d5_Y.loci
	// -purity       -pu   Purity (rho) setting for manual setting of sunrise plot location
    // -ploidy       -pi   Ploidy (psi) setting for manual setting of sunrise plot location
    // pu and pi are for plotting only, which we might not need. plus they have an auto setting. 
	species="HUMAN"
	assembly=37
	"""
	SPECIES="${species}" ; ASSEMBLY=${assembly} ; PROTOCOL="WGS"
	mkdir -p ascatResults ; export TMPDIR=\$(pwd)/tmp ; mkdir \$TMPDIR 
	ascat.pl \
	-o ./ascatResults \
	-t ${tumorBam} \
	-n ${normalBam} \
	-sg $REF_BASE/CNV_SV_ref/ascat/SnpGcCorrections.tsv \
	-r $REF_BASE/core_ref_GRCh37d5/genome.fa \
	-q 20 \
	-g L \
	-rs "\$SPECIES" \
	-ra \$ASSEMBLY \
	-pr \$PROTOCOL \
	-pl ILLUMINA \
	-c ${task.cpus} \
	-force 

	find . > listoffiles
	"""
}

process runBRASS {
	label 'BRASS'   
	input:
	set file(genomeFile), file(genomeIndex), file(REF_BASE) from Channel.value([file(genomeFile), file(genomeIndex), file(RefBase)])
	file("samplestatistics.txt") from ascatOut4Brass
	set file(tumorBam), file(tumorBai), file(normalBam), file(normalBai) from bamInput2
	set file(basFile1), file(basFile2) from basFiles

	output:
	file("listoffiles") into BRASSOutFiles

	script:
	species="HUMAN"
	assembly=37
	"""
	SPECIES="${species}" ; ASSEMBLY=${assembly} ; PROTOCOL="WGS"
	mkdir -p brassResults ; export TMPDIR=\$(pwd)/tmp ; mkdir \$TMPDIR 
	brass.pl -j 4 -k 4 -c ${task.cpus} \
	-d ${REF_BASE}/CNV_SV_ref/brass/HiDepth.bed.gz \
	-f ${REF_BASE}/CNV_SV_ref/brass/brass_np.groups.gz \
	-g ${REF_BASE}/core_ref_GRCh37d5/genome.fa \
	-s "\$SPECIES" -as \$ASSEMBLY -pr \$PROTOCOL -pl ILLUMINA \
	-g_cache ${REF_BASE}/VAGrENT_ref_GRCh37d5_ensembl_75/vagrent/vagrent.cache.gz \
	-vi ${REF_BASE}/CNV_SV_ref/brass/viral.genomic.fa.2bit \
	-mi ${REF_BASE}/CNV_SV_ref/brass/all_ncbi_bacteria \
	-b ${REF_BASE}/CNV_SV_ref/brass/500bp_windows.gc.bed.gz \
	-ct ${REF_BASE}/CNV_SV_ref/brass/CentTelo.tsv \
	-cb ${REF_BASE}/CNV_SV_ref/brass/cytoband.txt \
	-t $tumorBam \
	-n $normalBam \
	-ss samplestatistics.txt \
	-o brassResults
	find . > listoffiles
	"""
	
}

