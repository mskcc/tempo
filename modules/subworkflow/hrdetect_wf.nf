include { HRDetect }				from '../process/HRDetect/HRDetect'

workflow hrdetect_wf {
	take:
	CNVoutput
	MAFoutput 
	SVoutput
	
	main:
	HRDetectVariantsIn = MAFoutput
		.combine(CNVoutput,by:[0,1,2])
		.combine(SVoutput, by:[0,1,2])

	HRDetect(
		HRDetectVariantsIn,
		workflow.projectDir + "/containers/signaturetoolslib/HRDetect_wrapper.R"
	)
	
}
