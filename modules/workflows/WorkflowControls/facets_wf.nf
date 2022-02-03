include { DoFacets }                   from '../Facets/DoFacets' 
include { DoFacetsPreviewQC }          from '../Facets/DoFacetsPreviewQC' 

workflow facets_wf
{
  take:
    bamFiles

  main:
    referenceMap = params.referenceMap
    targetsMap   = params.targetsMap
    
    DoFacets(bamFiles, Channel.value([referenceMap.facetsVcf]))

    DoFacetsPreviewQC(DoFacets.out.Facets4FacetsPreview)

    DoFacets.out.FacetsRunSummary.combine(DoFacetsPreviewQC.out.FacetsPreviewOut, by:[0,1]).set{ FacetsQC4Aggregate }      // idTumor, idNormal, summaryFiles, qcFiles
    DoFacets.out.FacetsRunSummary.combine(DoFacetsPreviewQC.out.FacetsPreviewOut, by:[0,1]).set{ FacetsQC4SomaticMultiQC } // idTumor, idNormal, summaryFiles, qcFiles
    FacetsQC4Aggregate.map{ idTumor, idNormal, summaryFiles, qcFiles ->
      ["placeholder",idTumor, idNormal, summaryFiles, qcFiles]
    }.set{ FacetsQC4Aggregate }

  emit:
    snpPileupOutput         = DoFacets.out.snpPileupOutput
    FacetsOutput            = DoFacets.out.FacetsOutput
    FacetsOutLog4Aggregate  = DoFacets.out.FacetsOutLog4Aggregate
    FacetsPurity4Aggregate  = DoFacets.out.FacetsPurity4Aggregate
    FacetsHisens4Aggregate  = DoFacets.out.FacetsHisens4Aggregate
    facetsPurity            = DoFacets.out.facetsPurity
    facetsForMafAnno        = DoFacets.out.facetsForMafAnno
    Facets4FacetsPreview    = DoFacets.out.Facets4FacetsPreview
    FacetsArmGeneOutput     = DoFacets.out.FacetsArmGeneOutput
    FacetsArmLev4Aggregate  = DoFacets.out.FacetsArmLev4Aggregate
    FacetsGeneLev4Aggregate = DoFacets.out.FacetsGeneLev4Aggregate
    FacetsQC4MetaDataParser = DoFacets.out.FacetsQC4MetaDataParser
    FacetsRunSummary        = DoFacets.out.FacetsRunSummary
    FacetsPreviewOut        = DoFacetsPreviewQC.out.FacetsPreviewOut
    FacetsQC4Aggregate      = FacetsQC4Aggregate
    FacetsQC4SomaticMultiQC = FacetsQC4SomaticMultiQC
}