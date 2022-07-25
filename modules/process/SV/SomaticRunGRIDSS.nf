process SomaticRunGRIDSS {
	tag "${idTumor}__${idNormal}"
    publishDir "${params.outDir}/somatic/${idTumor}__${idNormal}/gridss", mode: params.publishDirMode, pattern: "*.{vcf.gz,vcf.gz.tbi}"

	input:
	tuple val(idTumor), val(idNormal), val(target), path(bamTumor), path(baiTumor), path(bamNormal), path(baiNormal)
	path(genomeFile)
    path(genomeIndex)
	path(bwaIndex)
	path(excludeRegions)
	path(breakend_pon)
	path(breakpoint_pon)
	path(breakpoint_hotspot)

	output:
	tuple val(idTumor), val(idNormal), val(target), path("${gripss_somatic_filtered_vcf}"), path("${gripss_somatic_filtered_vcf}.tbi")
	
	script:
	gridss_unfiltered_vcf="${idTumor}__${idNormal}.vcf.gz"
	gripss_somatic_vcf="${idTumor}__${idNormal}.somatic.vcf.gz"
	gripss_somatic_filtered_vcf="${idTumor}__${idNormal}.somatic.filtered.vcf.gz"
	gripss_jar="/opt/hmftools/gripss-1.9.jar"
	gridss_exec="/opt/gridss/gridss.sh"
	"""
	echo "chunkSize=50000000" > gridss.properties
	${gridss_exec} \\
	  --reference ${genomeFile} \\
	  --output ${idTumor}__${idNormal}.vcf.gz \\
	  --threads ${task.cpus * 2} \\
	  --labels ${idNormal},${idTumor} \\
	  --blacklist ${excludeRegions} \\
	  --assembly assembly.bam \\
	  --configuration gridss.properties \\
	  ${bamNormal} \\
	  ${bamTumor}

	java -Xmx24G -cp ${gripss_jar} com.hartwig.hmftools.gripss.GripssApplicationKt \\
	  -ref_genome ${genomeFile} \\
	  -breakpoint_hotspot ${breakpoint_hotspot} \\
	  -breakend_pon ${breakend_pon} \\
	  -breakpoint_pon ${breakpoint_pon} \\
	  -input_vcf ${gridss_unfiltered_vcf} \\
	  -output_vcf ${gripss_somatic_vcf} \\
	  -tumor ${idTumor} 2>&1 | tee ${idTumor}.gripss.log
	if [[ ! -f ${gripss_somatic_vcf} ]] ; then
	  echo "Error creating ${gripss_somatic_vcf}. Aborting"
	  exit 1
	fi
	
	java -Xmx24G -cp ${gripss_jar} com.hartwig.hmftools.gripss.GripssHardFilterApplicationKt \\
	  -input_vcf ${gripss_somatic_vcf} \\
	  -output_vcf ${gripss_somatic_filtered_vcf} 
	"""
}