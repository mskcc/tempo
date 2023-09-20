nextflow.enable.dsl=2
include { AggregateIndels }               from '../process/SNV/AggregateIndels' 


// include { defineReferenceMap; loadTargetReferences } from '../function/define_maps'

// include { touchInputs; watchMapping; watchBamMapping; watchPairing; watchAggregateWithResult; watchAggregate } from '../function/watch_inputs'
// include { defineReferenceMap; loadTargetReferences } from '../function/define_maps'

// include { touchInputs; watchMapping; watchBamMapping; watchPairing; watchAggregateWithResult; watchAggregate } from '../function/watch_inputs'

// pairingQc    = params.pairing
// referenceMap = defineReferenceMap()
// targetsMap   = loadTargetReferences()


workflow indel_wf {
  take: 
    bamFiles    
    strelkaOut
    platypusOut
    svabaIndelout

  main:
    // pairingQc    = params.pairing


    // refdir = Channel.value("/juno/work/ccs/orgeraj/indel_proj/pindelsetup/refarea/human_g1k_v37_decoy.fasta")
    // referenceMap = referenceMap

    refs=Channel.value(["/juno/work/ccs/orgeraj/indel_proj/pindelsetup/refarea/human_g1k_v37_decoy.fasta","/juno/work/taylorlab/cmopipeline/mskcc-igenomes/igenomes/Homo_sapiens/GATK/GRCh37/Annotation/GATKBundle/dbsnp_138.b37.vcf","/juno/work/ccs/orgeraj/indel_proj/pindelsetup/refarea/human_g1k_v37_decoy.fasta.fai"])

    combinedChannel = strelkaOut.combine(platypusOut, by:0)
    combinedChannelfin= combinedChannel.combine(svabaIndelout, by:0)
    combinedChannelfin.view()
   
    // refs.view()
    // filter to pass
    AggregateIndels(bamFiles,combinedChannelfin,refs)

    
    emit:
    tsvGroup  = AggregateIndels.out.tsvGroup
}

workflow {
    
    bamFiles = Channel.of(["P-0032509-T02-IM7__P-0032509-N01-IM7_impact505", "P-0032509-T02-IM7", "P-0032509-N01-IM7" , 
     "impact505", 
    "/juno/res/ci/share/ccs/delphi/chaos/delphiA/20230106_20_19_830974/bams/P-0032509-T02-IM7/P-0032509-T02-IM7.bam"	,
    "/juno/res/ci/share/ccs/delphi/chaos/delphiA/20230106_20_19_830974/bams/P-0032509-T02-IM7/P-0032509-T02-IM7.bam.bai",
    "/juno/res/ci/share/ccs/delphi/chaos/delphiA/20230106_20_19_830974/bams/P-0032509-N01-IM7/P-0032509-N01-IM7.bam",	
    "/juno/res/ci/share/ccs/delphi/chaos/delphiA/20230106_20_19_830974/bams/P-0032509-N01-IM7/P-0032509-N01-IM7.bam.bai"])

    strelkaOut=Channel.of(["P-0032509-T02-IM7__P-0032509-N01-IM7_impact505",
    "/juno/work/ccs/orgeraj/wf_indel_tempo/resultsImpact/somatic/P-0032509-T02-IM7__P-0032509-N01-IM7/strelka2/P-0032509-T02-IM7__P-0032509-N01-IM7.strelka2.vcf.gz",
    "/juno/work/ccs/orgeraj/wf_indel_tempo/resultsImpact/somatic/P-0032509-T02-IM7__P-0032509-N01-IM7/strelka2/P-0032509-T02-IM7__P-0032509-N01-IM7.strelka2.vcf.gz.tbi"])
    

    platypusOut=Channel.of(["P-0032509-T02-IM7__P-0032509-N01-IM7_impact505", 
    "/juno/work/ccs/orgeraj/wf_indel_tempo/resultsImpact/somatic/P-0032509-T02-IM7__P-0032509-N01-IM7/platypus/P-0032509-T02-IM7__P-0032509-N01-IM7.Somatic.Platypus.vcf"])
    
    svabaIndelout=Channel.of(["P-0032509-T02-IM7__P-0032509-N01-IM7_impact505", 
    "/juno/work/ccs/orgeraj/wf_indel_tempo/resultsImpact/somatic/P-0032509-T02-IM7__P-0032509-N01-IM7/svaba/P-0032509-T02-IM7__P-0032509-N01-IM7.reheader.svaba.somatic.indel.vcf.gz"])
    
    // refs=Channel.of(["/juno/work/ccs/orgeraj/indel_proj/pindelsetup/refarea/human_g1k_v37_decoy.fasta","/juno/work/taylorlab/cmopipeline/mskcc-igenomes/igenomes/Homo_sapiens/GATK/GRCh37/Annotation/GATKBundle/dbsnp_138.b37.vcf"])
    // referenceMap = defineReferenceMap()
    // targetsMap   = loadTargetReferences()
    // referenceMap = params.referenceMap

    indel_wf( bamFiles,strelkaOut,platypusOut,svabaIndelout )
}