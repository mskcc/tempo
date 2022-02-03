include { QcCollectHsMetrics }                 from '../QC/QcCollectHsMetrics' 
include { QcQualimap }                         from '../QC/QcQualimap' 
include { QcAlfred }                           from '../QC/QcAlfred'
include { SampleRunMultiQC }                   from '../QC/SampleRunMultiQC'
include { QcConpairAll }                       from '../QC/QcConpairAll'
include { CohortRunMultiQC }                   from '../QC/CohortRunMultiQC'

workflow sampleQC_wf
{
  take:
    inputChannel
    fastPJson

  main:
    referenceMap = params.referenceMap
    targetsMap   = params.targetsMap

    inputChannel.map{ idSample, target, bam, bai ->
        [idSample, target, bam, bai, targetsMap."$target".targetsInterval,  targetsMap."$target".baitsInterval]
    }.set{ bamsBQSR4HsMetrics }

    QcCollectHsMetrics(bamsBQSR4HsMetrics,
                       Channel.value([referenceMap.genomeFile, referenceMap.genomeIndex, referenceMap.genomeDict]))

    if (params.assayType != "exome"){
      QcCollectHsMetrics.out.collectHsMetricsOutput
        .map{ idSample, target, bam, bai, targetList, baitList -> [idSample, ""]}
        .set{ collectHsMetricsOutput }
    }

    inputChannel
      .map{ idSample, target, bam, bai -> [ idSample, target, bam, bai, file(targetsMap."$target".targetsBed) ]}
      .set{ bamsBQSR4Qualimap }

    QcQualimap(bamsBQSR4Qualimap)

    Channel.from(true, false).set{ ignore_read_groups }
    inputChannel
      .map{ idSample, target, bam, bai -> 
        [ idSample, target, bam, bai, targetsMap."$target".targetsBedGz, targetsMap."$target".targetsBedGzTbi ]
      }.set{ bamsBQSR4Alfred }

    QcAlfred(ignore_read_groups, 
             bamsBQSR4Alfred,
             Channel.value([referenceMap.genomeFile]))

    QcAlfred.out.alfredOutput
      .groupTuple(size:2, by:0)
      .join(fastPJson, by:0)
      .join(QcQualimap.out.qualimap4Process, by:0)
      .join(QcCollectHsMetrics.out.collectHsMetricsOutput, by:0)
      .set{ sampleMetrics4MultiQC }

    SampleRunMultiQC(sampleMetrics4MultiQC, 
                     Channel.value([params.multiqcWesConfig, params.multiqcWgsConfig, params.multiqcTempoLogo]))

  emit:
    bamsQcStats4Aggregate  = QcAlfred.out.bamsQcStats4Aggregate
    collectHsMetricsOutput = QcCollectHsMetrics.out.collectHsMetricsOutput
    qualimap4Process       = QcQualimap.out.qualimap4Process
}