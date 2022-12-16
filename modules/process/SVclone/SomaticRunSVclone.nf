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
  tuple val("placeHolder"), val(idTumor), val(idNormal),
    path("svclone/svs/*cluster_certainty.txt"),
    path("svclone/snvs/*cluster_certainty.txt"), emit: SVclone4Aggregate

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

mkdir -p svclone/svs svclone/snvs
cp ${outputPrefix}/ccube_out/post_assign/*.RData ${outputPrefix}/ccube_out/post_assign/*.pdf svclone/svs
cp ${outputPrefix}/ccube_out/post_assign/snvs/*.RData ${outputPrefix}/ccube_out/post_assign/snvs/*.pdf svclone/snvs
for i in ${outputPrefix}/ccube_out/post_assign/*.txt ; do
  sed "s/^/${outputPrefix}\\t/g" \$i | sed "0,/^${outputPrefix}\\t/s//sampleid\\t/" > svclone/svs/\$(basename \$i)
done
for i in ${outputPrefix}/ccube_out/post_assign/snvs/*.txt ; do
  sed "s/^/${outputPrefix}\\t/g" \$i | sed "0,/^${outputPrefix}\\t/s//sampleid\\t/" > svclone/snvs/\$(basename \$i)
done
"""
}
