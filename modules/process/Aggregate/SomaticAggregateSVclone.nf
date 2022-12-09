process SomaticAggregateSVclone {
	tag "${cohortid}"
    publishDir "${params.outDir}/cohort_level/${cohortid}", mode: params.publishDirMode

	input:
	tuple val(cohortid), path("sv/*"),path("snv/*")

	output:
	path("svclone_sv_cluster_certainty.tsv")
	path("svclone_snv_cluster_certainty.tsv")

	script:
	"""
	awk 'FNR==1 && NR!=1{next;}{print}' sv/* > svclone_sv_cluster_certainty.tsv
	awk 'FNR==1 && NR!=1{next;}{print}' snv/* > svclone_snv_cluster_certainty.tsv
	"""
}