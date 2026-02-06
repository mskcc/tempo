include { MetaDataParser }             from '../process/MetaParse/MetaDataParser' 

workflow mdParse_wf
{
  take:
    mergedChannelMetaDataParser

  main:
    MetaDataParser(mergedChannelMetaDataParser)

  emit:
    MetaData4Aggregate = MetaDataParser.out.MetaData4Aggregate
}
