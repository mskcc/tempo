process QcAlfred {
  tag {idSample + "@" + "ignore_rg_" + ignore_rg }

  publishDir "${params.outDir}/bams/${idSample}/alfred", mode: params.publishDirMode

  input:
    each ignore_rg
    tuple val(idSample), val(target), path(bam), path(bai), path(targets), path(targetsIndex)
    path(genomeFile)
    val(runQC)
    
  output:
    tuple val(idSample), path("${idSample}.alfred*tsv.gz"), emit: bamsQcStats4Aggregate
    tuple val(idSample), path("${idSample}.alfred*tsv.gz"), path("${idSample}.alfred*tsv.gz.pdf"), emit: alfredOutput

  when: runQC

  script:
  if (workflow.profile == "juno") {
    if (bam.size() > 200.GB) {
      task.time = { params.maxWallTime }
    }
    else if (bam.size() < 100.GB) {
      task.time = task.exitStatus.toString() in params.wallTimeExitCode ? { params.medWallTime } : { params.minWallTime }
    }
    else {
      task.time = task.exitStatus.toString() in params.wallTimeExitCode ? { params.maxWallTime } : { params.medWallTime }
    }
    task.time = task.attempt < 3 ? task.time : { params.maxWallTime }
  }

  options = ""
  if (params.assayType == "exome") {
    if (target == "agilent") options = "--bed ${targets}"
  }
  def ignore = ignore_rg ? "--ignore" : ""
  def outfile = ignore_rg ? "${idSample}.alfred.tsv.gz" : "${idSample}.alfred.per_readgroup.tsv.gz"
  """
  alfred qc ${options} \
    --reference ${genomeFile} \
    ${ignore} \
    --outfile ${outfile} \
    ${bam} && \
    Rscript --no-init-file /opt/alfred/scripts/stats.R ${outfile}
  """
}