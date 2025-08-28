include { MetaDataParser }             from '../process/MetaParse/MetaDataParser' 
include { defineReferenceMap; loadTargetReferences } from '../function/define_maps'

targetsMap   = loadTargetReferences()

workflow mdParse_wf
{
  take:
    inputSampleList
    facetsPurity
    maf4MetaDataParser
    FacetsQC4MetaDataParser
    msi4MetaDataParser
    mutSig4MetaDataParser
    hlaOutput

  main:
    mergedChannelMetaDataParser = inputSampleList
        .join(facetsPurity, by: [0,1,2], remainder: true)
        .join(maf4MetaDataParser, by: [0,1,2], remainder: true)
        .join(FacetsQC4MetaDataParser, by: [0,1,2], remainder: true)
        .join(msi4MetaDataParser, by: [0,1,2], remainder: true)
        .join(mutSig4MetaDataParser, by: [0,1,2], remainder: true)
        .join(hlaOutput, by: [1,2], remainder: true)
        .unique()
        .map { idNormal, target, idTumor, purityOut, mafFile, qcOutput, msifile, mutSig, placeHolder, polysolverFile ->
            def create_deterministic_placeholder = { purpose ->
                def source_placeholder = file(params.dummy_file)
                def unique_name = "${source_placeholder.baseName}_${purpose}.${source_placeholder.extension}"
                def unique_placeholder_copy = file(System.getProperty("java.io.tmpdir")).resolve(unique_name)
                if (!unique_placeholder_copy.exists()) {
                    unique_placeholder_copy.text = ''
                }
                return unique_placeholder_copy
            }

            def placeholder_val = ""

            def safe_purityOut = purityOut ?: create_deterministic_placeholder("purityOut")
            def safe_mafFile = mafFile ?: create_deterministic_placeholder("mafFile")
            def safe_qcOutput = qcOutput ?: create_deterministic_placeholder("qcOutput")
            def safe_msifile = msifile ?: create_deterministic_placeholder("msifile")
            def safe_mutSig = mutSig ?: create_deterministic_placeholder("mutSig")
            def safe_polysolverFile = polysolverFile ?: create_deterministic_placeholder("polysolverFile")

            def safe_placeHolder = placeHolder ?: placeholder_val

            def codingBed = targetsMap[target]?.codingBed ? file(targetsMap[target].codingBed) : create_deterministic_placeholder("codingBed")

            return tuple(
                idNormal,
                target,
                idTumor,
                safe_purityOut,
                safe_mafFile,
                safe_qcOutput,
                safe_msifile,
                safe_mutSig,
                safe_placeHolder,
                safe_polysolverFile,
                codingBed
            )
        }
    MetaDataParser(mergedChannelMetaDataParser)

  emit:
    MetaData4Aggregate = MetaDataParser.out.MetaData4Aggregate
}
