process QcQualimap {
  tag {idSample}
  
  publishDir "${params.outDir}/bams/${idSample}/qualimap", mode: params.publishDirMode, pattern: "*.{html,tar.gz}"
  publishDir "${params.outDir}/bams/${idSample}/qualimap", mode: params.publishDirMode, pattern: "*/*"

  input:
    tuple val(idSample), val(target), path(bam), path(bai), path(targetsBed)
    val(runQC)

  output:
    tuple val(idSample), path("${idSample}_qualimap_rawdata.tar.gz"), emit: qualimap4Process
    tuple val(idSample), path("*.html"), path("css/*"), path("images_qualimapReport/*"), emit: qualimapOutput
  
  when: runQC   

  script:
  if (params.assayType == "exome"){
    gffOptions = "-gff ${targetsBed}"
    nr = 750
    nw = 300
  } else { 
    gffOptions = "-gd HUMAN" 
    nr = 500
    nw = 300
  }
  availMem = task.cpus * task.memory.toString().split(" ")[0].toInteger()
  // javaMem = availMem > 20 ? availMem - 4 : ( availMem > 10 ? availMem - 2 : ( availMem > 1 ? availMem - 1 : 1 ))
  javaMem = availMem > 20 ? (availMem * 0.75).round() : ( availMem > 1 ? availMem - 1 : 1 )
  if (workflow.profile == "juno") {
    if (bam.size() > 200.GB) {
      task.time = { params.maxWallTime }
    }
    else if (bam.size() < 100.GB) {
      task.time = task.exitStatus.toString() in params.wallTimeExitCode ? { params.medWallTime } : { params.minWallTime }
    }
    else {
      task.time = task.exitStatus.toString() in params.wallTimeExitCode ? { params.maxWallTime } : { params.medWallTime }
    }
    task.time = task.attempt < 3 ? task.time : { params.maxWallTime }
  }
  """
  qualimap bamqc \
  -bam ${bam} \
  ${gffOptions} \
  -outdir ${idSample} \
  -nt ${ task.cpus * 2 } \
  -nw ${nw} \
  -nr ${nr} \
  --java-mem-size=${javaMem}G

  mv ${idSample}/* . 
  tar -czf ${idSample}_qualimap_rawdata.tar.gz genome_results.txt raw_data_qualimapReport/* 
  """
}