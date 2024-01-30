include { DoFacets }                   from '../process/Facets/DoFacets' 
include { DoFacetsPreviewQC }          from '../process/Facets/DoFacetsPreviewQC' 

workflow facets_wf
{
  take:
    bamFiles

  main:
    referenceMap = params.referenceMap
    targetsMap   = params.targetsMap
    outputDir = "facets${params.facets.R_lib}c${params.facets.cval}pc${params.facets.purity_cval}"
    
    DoFacets(
      bamFiles, 
      referenceMap.facetsVcf,
      workflow.projectDir + "/containers/facets-suite-preview-htstools",  
      outputDir
    )

    DoFacetsPreviewQC(DoFacets.out.Facets4FacetsPreview)

    DoFacets.out.FacetsRunSummary.combine(DoFacetsPreviewQC.out.FacetsPreviewOut, by:[0,1]).set{ FacetsQC4Aggregate }      // idTumor, idNormal, summaryFiles, qcFiles
    DoFacets.out.FacetsRunSummary.combine(DoFacetsPreviewQC.out.FacetsPreviewOut, by:[0,1]).set{ FacetsQC4SomaticMultiQC } // idTumor, idNormal, summaryFiles, qcFiles
    FacetsQC4Aggregate.map{ idTumor, idNormal, summaryFiles, qcFiles ->
      ["placeholder",idTumor, idNormal, summaryFiles, qcFiles]
    }.set{ FacetsQC4Aggregate }

    DoFacets.out.facets4Aggregate
      .map{
        ["placeHolder"] + it
      }.set{facets4Aggregate}

    DoFacets.out.FacetsArmGeneOutput
      .map{
        ["placeHolder"] + it
      }.set{FacetsArmGeneOutput}

  emit:
    snpPileupOutput            = DoFacets.out.snpPileupOutput
    FacetsOutput               = DoFacets.out.FacetsOutput
    facets4Aggregate           = facets4Aggregate
    facetsPurity               = DoFacets.out.facetsPurity
    facetsForMafAnno           = DoFacets.out.facetsForMafAnno
    Facets4FacetsPreview       = DoFacets.out.Facets4FacetsPreview
    FacetsArmGeneOutput        = FacetsArmGeneOutput
    FacetsQC4MetaDataParser    = DoFacets.out.FacetsQC4MetaDataParser
    FacetsRunSummary           = DoFacets.out.FacetsRunSummary
    FacetsPreviewOut           = DoFacetsPreviewQC.out.FacetsPreviewOut
    FacetsQC4Aggregate         = FacetsQC4Aggregate
    FacetsQC4SomaticMultiQC    = FacetsQC4SomaticMultiQC
    FacetsHisensCNV4HrDetect            = DoFacets.out.FacetsHisensCNV4HrDetect
    FacetsHisensCNV4HrDetectFiltered    = DoFacets.out.FacetsHisensCNV4HrDetectFiltered
    FacetsHisensSampleStatistics4BRASS  = DoFacets.out.FacetsHisensSampleStatistics4BRASS
    FacetsPurityCNV4HrDetect            = DoFacets.out.FacetsPurityCNV4HrDetect
    FacetsPurityCNV4HrDetectFiltered    = DoFacets.out.FacetsPurityCNV4HrDetectFiltered
    FacetsPuritySampleStatistics4BRASS  = DoFacets.out.FacetsPuritySampleStatistics4BRASS

}
