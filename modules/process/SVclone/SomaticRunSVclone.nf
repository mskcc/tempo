process SomaticRunSVclone {
tag "${idTumor}__${idNormal}"

publishDir "${params.outDir}/somatic/${outputPrefix}/clonality/", mode: params.publishDirMode

input:
  tuple val(idTumor), val(idNormal), val(target), 
    path(bamTumor), path(baiTumor),
    path(bamNormal), path(baiNormal),
    path(inBedpe),
    path(mafFiltered),
    path(cnv),
    path(ploidyIn)
  path(prepare_svclone_input_script)
output:
  tuple val(idTumor), val(idNormal), val(target),
    path("${outputPrefix}")

script:
outputPrefix = "${idTumor}__${idNormal}"
"""
echo "Preparing svclone inputs using custom script"
python ${prepare_svclone_input_script} \\
  --cfg_template /config/svclone_config.ini \\
  --bedpe ${inBedpe} \\
  --maf ${mafFiltered} \\
  --purity_ploidy ${ploidyIn} \\
  --out_dir svclone_in \\
  --sampleid ${outputPrefix} \\
  --bam ${bamTumor}

echo "Running SVclone"
echo "Running svclone annotate"
svclone annotate \\
  -i svclone_in/simple.sv.txt \\
  -b ${bamTumor} \\
  -s ${outputPrefix} \\
  --sv_format simple \\
  -cfg svclone_in/svclone_config.ini

echo "Running svclone count"
svclone count \\
  -i ${outputPrefix}/${outputPrefix}_svin.txt \\
  -b ${bamTumor} \\
  -cfg svclone_in/svclone_config.ini \\
  -s ${outputPrefix}

echo "Running svclone filter"
svclone filter \\
  -s ${outputPrefix} \\
  -i ${outputPrefix}/${outputPrefix}_svinfo.txt \\
  -p svclone_in/svclone_ploidy.txt \\
  -cfg svclone_in/svclone_config.ini \\
  --cnvs ${cnv} \\
  --snvs svclone_in/callstats.txt \\
  --snv_format mutect_callstats

echo "Running svclone cluster"
svclone cluster \\
  -cfg svclone_in/svclone_config.ini \\
  -s ${outputPrefix}

echo "Running svclone postassign"
svclone postassign \\
  -s ${outputPrefix} \\
  --joint

"""
}
