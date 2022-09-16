include { SomaticRunSVclone }          from '../process/SVclone/SomaticRunSVclone'

workflow clonality_wf {
	take:
	bams
	final_bedpe
	final_maf
	final_cnv
	final_purity_ploidy
	
	main:
	SVcloneInput = bams
		.combine(final_bedpe,by:[0,1,2])
		.combine(final_maf, by:[0,1,2])
		.combine(final_cnv, by:[0,1,2])
		.combine(final_purity_ploidy, by:[0,1,2])

	SomaticRunSVclone(
		SVcloneInput,
		workflow.projectDir + "/containers/svclone/prepare_svclone_inputs.py"
	)
	
}
