process RunNeoantigen {
  tag "${idTumor + "__" + idNormal}"

  publishDir "${params.outDir}/somatic/${outputPrefix}/neoantigen/", mode: params.publishDirMode, pattern: "*.txt"

  input:
    tuple val(idNormal), val(target), val(placeHolder), path(polysolverFile), val(idTumor), path(mafFile)
    tuple path(neoantigenCDNA), path(neoantigenCDS)

  output:
    tuple val(placeHolder), val(idTumor), val(idNormal), path("${idTumor}__${idNormal}.all_neoantigen_predictions.txt"), emit: NetMhcStats4Aggregate
    path("${idTumor}__${idNormal}.all_neoantigen_predictions.txt"), emit: NetMhcStatsOutput
    tuple val(idTumor), val(idNormal), val(target), path("${outputDir}/${outputPrefix}.neoantigens.maf"), emit: mafFileForMafAnno

  script:

  if (workflow.profile == "juno") {
    if(mafFile.size() > 10.MB){
      task.time = { params.maxWallTime }
    }
    else if (mafFile.size() < 5.MB){
      task.time = task.exitStatus.toString() in params.wallTimeExitCode ? { params.medWallTime } : { params.minWallTime }
    }
    else {
      task.time = task.exitStatus.toString() in params.wallTimeExitCode ? { params.maxWallTime } : { params.medWallTime }
    }
    task.time = task.attempt < 3 ? task.time : { params.maxWallTime }
  }

  outputPrefix = "${idTumor}__${idNormal}"
  outputDir = "neoantigen"
  tmpDir = "${outputDir}-tmp"
  tmpDirFullPath = "\$PWD/${tmpDir}/"  // must set full path to tmp directories for netMHC and netMHCpan to work; for some reason doesn't work with /scratch, so putting them in the process workspace
  """
  export TMPDIR=${tmpDirFullPath}
  mkdir -p ${tmpDir}
  chmod 777 ${tmpDir}

  python /usr/local/bin/neoantigen/neoantigen.py \
    --config_file /usr/local/bin/neoantigen/neoantigen-docker.config \
    --sample_id ${outputPrefix} \
    --hla_file ${polysolverFile} \
    --maf_file ${mafFile} \
    --output_dir ${outputDir}

  awk 'NR==1 {printf("%s\\t%s\\n", "sample", \$0)} NR>1 {printf("%s\\t%s\\n", "${outputPrefix}", \$0) }' neoantigen/*.all_neoantigen_predictions.txt > ${outputPrefix}.all_neoantigen_predictions.txt
  """
}
