process SomaticRunSVclone {
tag "${idTumor}__${idNormal}"

publishDir "${params.outDir}/somatic/${outputPrefix}/", mode: params.publishDirMode, pattern: "svclone/*"

input:
  tuple val(idTumor), val(idNormal), val(target), 
    path(bamTumor), path(baiTumor),
    path(bamNormal), path(baiNormal),
    path(inBedpe),
    path(mafFiltered),
    path(cnv),
    path(ploidyIn)
  path(svclone_wrapper)
output:
  tuple val(idTumor), val(idNormal), val(target),
    path("${outputPrefix}"), emit: SVcloneOutput
  tuple val(idTumor), val(idNormal), val(target),
    path("svclone/*"), emit: SVclonePublish

script:
outputPrefix = "${idTumor}__${idNormal}"
"""
python ${svclone_wrapper} \\
  --cfg_template /config/svclone_config.ini \\
  --bedpe ${inBedpe} \\
  --maf ${mafFiltered} \\
  --purity_ploidy ${ploidyIn} \\
  --out_dir svclone_in \\
  --sampleid ${outputPrefix} \\
  --bam ${bamTumor} \\
  --cnv ${cnv}

mkdir -p svclone/svs
cp -R ${outputPrefix}/ccube_out/post_assign/* svclone
mv svclone/*.txt svclone/*.RData svclone/*.pdf svclone/svs
"""
}
