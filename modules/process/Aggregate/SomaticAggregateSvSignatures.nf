process SomaticAggregateSvSignatures {
  tag { cohortID }
  publishDir "${params.outDir}/cohort_level/${cohortID}", mode: params.publishDirMode

input:
  tuple val(cohortID),
    path(svCataloguePdf),
    path(exposureFiles)

output:
  path("sv_catalogues.pdf") 
  path("sv_exposures.tsv")

script:
  """
  awk 'FNR==1 && NR!=1{next;}{print}' *_exposures.tsv > sv_exposures.tsv
  gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -dPDFSETTINGS=/prepress -sOutputFile=new.pdf *_catalogues.pdf
  mv new.pdf sv_catalogues.pdf
  """
}
