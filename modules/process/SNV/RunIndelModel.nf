process RunIndelModel {
  tag "${idTumor + "__" + idNormal}"

  publishDir "${params.outDir}/somatic/${outputPrefix}/IndelModelWF", mode: params.publishDirMode

  input:
    tuple val(idTumor), val(idNormal), val(target), path(bamTumor), path(baiTumor), path(bamNormal), path(baiNormal), path(inputTSV)
    tuple path(genomeFile), path(genomeIndex), path(repeatMaskerFile), path(py2bitfile)

  output:
    tuple val(idTumor), val(idNormal), val(target), path("${outputPrefix}_model.tsv"),  path("${outputPrefix}_IndelModelPass.vcf"), emit: indelOut 

  script:
  outputPrefix = "${idTumor}__${idNormal}"

  """

    gunzip -c ${repeatMaskerFile} > rmsk_mod.bed


    python3 /usr/bin/Eval_Indels.py \
    --indel-file ${inputTSV} \
    --rmsk rmsk_mod.bed \
    --model /Models/AUC_loss_1LSTM_simple_1206_99e_compiled_excellent/ \
    --py2bit ./${py2bitfile} \
    --outfile ${outputPrefix}_model.tsv


    python3 /usr/bin/Modified_tsv2vcf.py \
    -tsv ${outputPrefix}_model.tsv \
    -vcf ${outputPrefix}_IndelModelPass.vcf \
    -pass .8 \
    -tools Strelka Platypus Svaba \
    -paired  

    """
}