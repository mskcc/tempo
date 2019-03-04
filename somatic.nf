#!/usr/bin/env nextflow

/*
================================================================================
--------------------------------------------------------------------------------
 Processes overview
 - dellyCall
 - dellyFilter
*/


/*
================================================================================
=                           C O N F I G U R A T I O N                          =
================================================================================
*/

tsvPath = ''
if (params.sample) tsvPath = params.sample

referenceMap = defineReferenceMap()

bamFiles = Channel.empty()

tsvFile = file(tsvPath)

bamFiles = extractBamFiles(tsvFile)

( bamsForDelly, bamsForMutect2, bamsForManta, bamsForStrelka ) = bamFiles.into(4)

/*
================================================================================
=                               P R O C E S S E S                              =
================================================================================
*/

process dellyCall {
  tag { "DELLYCALL_" + idTumor + "_" + idNormal }

  publishDir "${ params.outDir }/VariantCalling/delly_call"

  input:
    set idTumor, idNormal, file(bamTumor), file(bamNormal), file(baiTumor), file(baiNormal) from bamsForDelly
    set file(genomeFile), file(genomeIndex), file(genomeDict), file(dbsnp), file(dbsnpIndex), file(knownIndels), file(knownIndelsIndex), file(intervals) from Channel.value([
      referenceMap.genomeFile,
      referenceMap.genomeIndex,
      referenceMap.genomeDict,
      referenceMap.dbsnp,
      referenceMap.dbsnpIndex,
      referenceMap.knownIndels,
      referenceMap.knownIndelsIndex,
      referenceMap.intervals,
    ])

  output:
    set file("${idTumor}_${idNormal}_DUP.bcf"), file("${idTumor}_${idNormal}_DUP.bcf.csi"),
        file("${idTumor}_${idNormal}_DEL.bcf"), file("${idTumor}_${idNormal}_DEL.bcf.csi"),
        file("${idTumor}_${idNormal}_INS.bcf"), file("${idTumor}_${idNormal}_INS.bcf.csi"),
        file("${idTumor}_${idNormal}_INV.bcf"), file("${idTumor}_${idNormal}_INV.bcf.csi"),
        file("${idTumor}_${idNormal}_BND.bcf"), file("${idTumor}_${idNormal}_BND.bcf.csi") into dellyCallOutput
    set idTumor, idNormal into forMakingSampleFile

  """
  sv_variants=("DUP" "BND" "DEL" "INS" "INV") 
  for sv_variant in "\${sv_variants[@]}"; 
  do
    outfile="${idTumor}_${idNormal}_\${sv_variant}.bcf" 
    delly call \
      -t "\${sv_variant}" \
      -o "\${outfile}" \
      -g ${genomeFile} \
      ${bamTumor} \
      ${bamNormal}
  done
  """
}

process makeSamplesFile {
  tag { "SAMPLESFILE_" + idTumor + "_" + idNormal }

  input:
    set idTumor, idNormal from forMakingSampleFile 

  output:
    set idTumor, idNormal, file("samples.tsv") into sampleTSVFile

  """
    echo "${idTumor}\ttumor\n${idNormal}\tcontrol" > samples.tsv
  """
} 

process dellyFilter {
  tag { "DELLYFILTER_" + idTumor + "_" + idNormal }

  publishDir "${ params.outDir }/VariantCalling/delly_filter"

  input:
    set idTumor, idNormal, file("samples.tsv") from sampleTSVFile 
    set file("${idTumor}_${idNormal}_DUP.bcf"), file("${idTumor}_${idNormal}_DUP.bcf.csi"),
        file("${idTumor}_${idNormal}_DEL.bcf"), file("${idTumor}_${idNormal}_DEL.bcf.csi"),
        file("${idTumor}_${idNormal}_INS.bcf"), file("${idTumor}_${idNormal}_INS.bcf.csi"),
        file("${idTumor}_${idNormal}_INV.bcf"), file("${idTumor}_${idNormal}_INV.bcf.csi"),
        file("${idTumor}_${idNormal}_BND.bcf"), file("${idTumor}_${idNormal}_BND.bcf.csi") from dellyCallOutput

  output:
    set file("${idTumor}_${idNormal}_DUP.filter.bcf"), file("${idTumor}_${idNormal}_DUP.filter.bcf.csi"),
        file("${idTumor}_${idNormal}_DEL.filter.bcf"), file("${idTumor}_${idNormal}_DEL.filter.bcf.csi"),
        file("${idTumor}_${idNormal}_INS.filter.bcf"), file("${idTumor}_${idNormal}_INS.filter.bcf.csi"),
        file("${idTumor}_${idNormal}_INV.filter.bcf"), file("${idTumor}_${idNormal}_INV.filter.bcf.csi"),
        file("${idTumor}_${idNormal}_BND.filter.bcf"), file("${idTumor}_${idNormal}_BND.filter.bcf.csi") into dellyFilterOutput 

  """
  sv_variants=("DUP" "BND" "DEL" "INS" "INV") 
  for sv_variant in "\${sv_variants[@]}"; 
  do
    delly_call_file="${idTumor}_${idNormal}_\${sv_variant}.bcf" 
    outfile="${idTumor}_${idNormal}_\${sv_variant}.filter.bcf" 
    delly filter \
      -f somatic \
      -o "\${outfile}" \
      -s "samples.tsv" \
      "\${delly_call_file}"
  done
  """
}

process runMutect2 {
  tag {"MUTECT2_" + idTumor + "_" + idNormal }

  publishDir "${ params.outDir }/VariantCalling/mutect2"

  input:
    set idTumor, idNormal, file(bamTumor), file(bamNormal), file(baiTumor), file(baiNormal) from bamsForMutect2
    set file(genomeFile), file(genomeIndex), file(genomeDict), file(dbsnp), file(dbsnpIndex) from Channel.value([
      referenceMap.genomeFile,
      referenceMap.genomeIndex,
      referenceMap.genomeDict,
      referenceMap.dbsnp,
      referenceMap.dbsnpIndex
    ])

  output:
    set idNormal, idTumor, file("${idTumor}_vs_${idNormal}_somatic.vcf") into mutect2Output

  script:
  """
  # Xmx hard-coded for now due to lsf bug
  gatk --java-options "-Xmx8g" \
    Mutect2 \
    -R ${genomeFile}\
    -I ${bamTumor}  -tumor ${idTumor} \
    -I ${bamNormal} -normal ${idNormal} \
    -O "${idTumor}_vs_${idNormal}_somatic.vcf"
  """
}

process runManta {
  tag {"RUNMANTA_" + idTumor + "_" + idNormal}

  publishDir "${params.outDir}/VariantCalling/Manta"

  input:
    set idTumor, idNormal, file(bamTumor), file(bamNormal), file(baiTumor), file(baiNormal) from bamsForManta
    set file(genomeFile), file(genomeIndex) from Channel.value([
      referenceMap.genomeFile,
      referenceMap.genomeIndex
    ])

  output:
    set idNormal, idTumor, file("*.vcf.gz"), file("*.vcf.gz.tbi") into mantaOutput
    set idNormal, idTumor, file("*.candidateSmallIndels.vcf.gz"), file("*.candidateSmallIndels.vcf.gz.tbi") into mantaToStrelka

  script:
  """
  configManta.py \
  --normalBam ${bamNormal} \
  --tumorBam ${bamTumor} \
  --reference ${genomeFile} \
  --runDir Manta

  python Manta/runWorkflow.py -m local -j ${task.cpus}

  mv Manta/results/variants/candidateSmallIndels.vcf.gz \
    Manta_${idTumor}_vs_${idNormal}.candidateSmallIndels.vcf.gz
  mv Manta/results/variants/candidateSmallIndels.vcf.gz.tbi \
    Manta_${idTumor}_vs_${idNormal}.candidateSmallIndels.vcf.gz.tbi
  mv Manta/results/variants/candidateSV.vcf.gz \
    Manta_${idTumor}_vs_${idNormal}.candidateSV.vcf.gz
  mv Manta/results/variants/candidateSV.vcf.gz.tbi \
    Manta_${idTumor}_vs_${idNormal}.candidateSV.vcf.gz.tbi
  mv Manta/results/variants/diploidSV.vcf.gz \
    Manta_${idTumor}_vs_${idNormal}.diploidSV.vcf.gz
  mv Manta/results/variants/diploidSV.vcf.gz.tbi \
    Manta_${idTumor}_vs_${idNormal}.diploidSV.vcf.gz.tbi
  mv Manta/results/variants/somaticSV.vcf.gz \
    Manta_${idTumor}_vs_${idNormal}.somaticSV.vcf.gz
  mv Manta/results/variants/somaticSV.vcf.gz.tbi \
    Manta_${idTumor}_vs_${idNormal}.somaticSV.vcf.gz.tbi
  """
}

// Running Strelka Best Practice with Manta indel candidates
// For easier joining, remaping channels to idTumor, idNormal...

bamsForStrelka = bamsForStrelka.map {
  idTumor, idNormal, bamTumor, bamNormal, baiTumor, baiNormal ->
  [idTumor, idNormal, bamTumor, bamNormal, baiTumor, baiNormal]
}.join(mantaToStrelka, by:[0,1,2]).map {
  idTumor, idNormal, bamTumor, bamNormal, baiTumor, baiNormal, mantaCSI, mantaCSIi ->
  [idTumor, idNormal, bamTumor, bamNormal, baiTumor, baiNormal, mantaCSI, mantaCSIi]
}

process runStrelka {
  tag {"RUNSTRELKA_" + idTumor + "_" + idNormal}

  publishDir "${params.outDir}/VariantCalling/Strelka"

  input:
    set  idTumor, idNormal, file(bamTumor), file(bamNormal), file(baiTumor), file(baiNormal), file(mantaCSI), file(mantaCSIi) from bamsForStrelka
    set file(genomeFile), file(genomeIndex), file(genomeDict) from Channel.value([
      referenceMap.genomeFile,
      referenceMap.genomeIndex,
      referenceMap.genomeDict
    ])

  output:
    set idNormal, idTumor, file("*.vcf.gz"), file("*.vcf.gz.tbi") into strelkaOutput


  script:
  """
  configureStrelkaSomaticWorkflow.py \
  --tumor ${bamTumor} \
  --normal ${bamNormal} \
  --referenceFasta ${genomeFile} \
  --indelCandidates ${mantaCSI} \
  --runDir Strelka

  python Strelka/runWorkflow.py -m local -j ${task.cpus}

  mv Strelka/results/variants/somatic.indels.vcf.gz \
    Strelka_${idTumor}_vs_${idNormal}_somatic_indels.vcf.gz
  mv Strelka/results/variants/somatic.indels.vcf.gz.tbi \
    Strelka_${idTumor}_vs_${idNormal}_somatic_indels.vcf.gz.tbi
  mv Strelka/results/variants/somatic.snvs.vcf.gz \
    Strelka_${idTumor}_vs_${idNormal}_somatic_snvs.vcf.gz
  mv Strelka/results/variants/somatic.snvs.vcf.gz.tbi \
    Strelka_${idTumor}_vs_${idNormal}_somatic_snvs.vcf.gz.tbi
  """
}

/*
================================================================================
=                               AWESOME FUNCTIONS                             =
================================================================================
*/

def checkParamReturnFile(item) {
  params."${item}" = params.genomes[params.genome]."${item}"
  return file(params."${item}")
}

def defineReferenceMap() {
  if (!(params.genome in params.genomes)) exit 1, "Genome ${params.genome} not found in configuration"
  return [
    'dbsnp'            : checkParamReturnFile("dbsnp"),
    'dbsnpIndex'       : checkParamReturnFile("dbsnpIndex"),
    // genome reference dictionary
    'genomeDict'       : checkParamReturnFile("genomeDict"),
    // FASTA genome reference
    'genomeFile'       : checkParamReturnFile("genomeFile"),
    // genome .fai file
    'genomeIndex'      : checkParamReturnFile("genomeIndex"),
    // BWA index files
    'bwaIndex'         : checkParamReturnFile("bwaIndex"),
    // intervals file for spread-and-gather processes
    'intervals'        : checkParamReturnFile("intervals"),
    // VCFs with known indels (such as 1000 Genomes, Mill’s gold standard)
    'knownIndels'      : checkParamReturnFile("knownIndels"),
    'knownIndelsIndex' : checkParamReturnFile("knownIndelsIndex"),
    // for SNP Pileup
    'facetsVcf'        : checkParamReturnFile("facetsVcf"),
    'facetsVcfIndex'   : checkParamReturnFile("facetsVcfIndex")
  ]
}

def extractBamFiles(tsvFile) {
  // Channeling the TSV file containing FASTQ.
  // Format is: "idTumor idNormal bamTumor bamNormal baiTumor baiNormal"
  Channel.from(tsvFile)
  .splitCsv(sep: '\t')
  .map { row ->
    SarekUtils.checkNumberOfItem(row, 6)
    def idTumor = row[0]
    def idNormal = row[1]
    def bamTumor = SarekUtils.returnFile(row[2])
    def bamNormal = SarekUtils.returnFile(row[3])
    def baiTumor = SarekUtils.returnFile(row[4])
    def baiNormal = SarekUtils.returnFile(row[5])

    SarekUtils.checkFileExtension(bamTumor,".bam")
    SarekUtils.checkFileExtension(bamNormal,".bam")
    SarekUtils.checkFileExtension(baiTumor,".bai")
    SarekUtils.checkFileExtension(baiNormal,".bai")

    [ idTumor, idNormal, bamTumor, bamNormal, baiTumor, baiNormal ]
  }
}
