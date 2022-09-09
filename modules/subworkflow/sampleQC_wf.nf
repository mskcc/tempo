include { QcCollectHsMetrics }                 from '../process/QC/QcCollectHsMetrics' 
include { QcQualimap }                         from '../process/QC/QcQualimap' 
include { QcAlfred }                           from '../process/QC/QcAlfred'
include { SampleRunMultiQC }                   from '../process/QC/SampleRunMultiQC'
include { QcConpairAll }                       from '../process/QC/QcConpairAll'

workflow sampleQC_wf
{
  take:
    inputChannel
    fastPJson
    intervals
    targetsBed
    targetsBedGz

  main:
    referenceMap = params.referenceMap
    targetsMap   = params.targetsMap

    inputChannel.combine(intervals)
	.filter{ idSample, target, bam, bai, target2, targetsInterval, baitsInterval ->
	    target == target2
	}.map{ idSample, target, bam, bai, target2, targetsInterval, baitsInterval ->
	    [idSample, target, bam, bai, targetsInterval, baitsInterval]
	}.set{ bamsBQSR4HsMetrics }

    QcCollectHsMetrics(bamsBQSR4HsMetrics,
                       Channel.value([referenceMap.genomeFile, referenceMap.genomeIndex, referenceMap.genomeDict]))

    if (params.assayType == "exome"){
	collectHsMetricsOutput = QcCollectHsMetrics.out.collectHsMetricsOutput
    } else {
    	inputChannel
        	.map{ idSample, target, bam, bai -> [idSample, ""]}
        	.set{ collectHsMetricsOutput }
    }

    inputChannel
      .combine(targetsBed)
      .filter{ idSample, target, bam, bai, targets2, bedfile -> target == targets2 }
      .map{ idSample, target, bam, bai, targets2, bedfile -> [idSample, target, bam, bai, bedfile] }
      .set{ bamsBQSR4Qualimap }

    QcQualimap(bamsBQSR4Qualimap)

    Channel.from(true, false).set{ ignore_read_groups }
    inputChannel
      .combine(targetsBedGz)
      .filter{ idSample, target, bam, bai, targets2, bedGz, bedGzTbi -> target == targets2 }
      .map{ idSample, target, bam, bai, targets2, bedGz, bedGzTbi -> [idSample, target, bam, bai, bedGz, bedGzTbi] }
      .set{ bamsBQSR4Alfred }

    QcAlfred(ignore_read_groups, 
             bamsBQSR4Alfred,
             Channel.value([referenceMap.genomeFile]))

    QcAlfred.out.alfredOutput
      .groupTuple(size:2, by:0)
      .join(fastPJson, by:0)
      .join(QcQualimap.out.qualimap4Process, by:0)
      .join(collectHsMetricsOutput, by:0)
      .set{ sampleMetrics4MultiQC }

    SampleRunMultiQC(sampleMetrics4MultiQC, 
                     Channel.value([params.multiqcWesConfig, params.multiqcWgsConfig, params.multiqcTempoLogo]))

  emit:
    bamsQcStats4Aggregate  = QcAlfred.out.bamsQcStats4Aggregate
    collectHsMetricsOutput = collectHsMetricsOutput 
    qualimap4Process       = QcQualimap.out.qualimap4Process
}
