include { GermlineFacetsAnnotation }          from '../GermSNV/GermlineFacetsAnnotation' 

workflow germlineSNV_facets
{
  take:
    facetsForMafAnno
    mafFileGermline

  main:
    facetsForMafAnno.combine(mafFileGermline, by: [0,1,2])
        .set{ facetsMafFileGermline }

    GermlineFacetsAnnotation(facetsMafFileGermline)

  emit:
    mafFile4AggregateGermline = GermlineFacetsAnnotation.out.mafFile4AggregateGermline
}