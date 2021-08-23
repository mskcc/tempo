process CreateScatteredIntervals {
  tag {targetId}

  input:
    tuple path(genomeFile), path(genomeIndex), path(genomeDict)
    tuple val(targetId), path(targets), path(targetsIndex)
    val(runSomatic)
    val(runGermline)
    
  output:
    tuple path("*.interval_list"), val(targetId), val(targetId), emit: mergedIList

  when: runSomatic || runGermline

  script:
  scatterCount = params.scatterCount
  """
  gatk SplitIntervals \
    --reference ${genomeFile} \
    --intervals ${targets} \
    --scatter-count ${scatterCount} \
    --subdivision-mode BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW \
    --output $targetId

  for i in $targetId/*.interval_list;
  do
    BASENAME=`basename \$i`
    mv \$i ${targetId}-\$BASENAME
  done
  """
}
