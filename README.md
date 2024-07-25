[![Build Status](https://travis-ci.com/mskcc/tempo.svg?token=DokCkCiDp43sqzeuXUHD&branch=master)](https://travis-ci.com/mskcc/tempo)

# <img id="logo" src="./docs/tempoLogo.png" height="42" align="left"/> TEMPO

Tempo is a CMO Computational Sciences (CCS) research pipeline processing WES & WGS tumor-normal pairs using the [Nextflow framework](https://www.nextflow.io/). Currently the pipeline is composed of alignment and QC, and detection of both somatic alterations and germline variants. Users can begin with inputs of either paired-end FASTQs or BAMs, and process these via the command line. 

For further details of how to begin processing data with Tempo, please view our [documentation](https://cmotempo.netlify.com/). For contributing to this project, please make a pull request as detailed [here](https://cmotempo.netlify.com/contributing-to-tempo.html).

The inspiration for this project derives from [Sarek](https://github.com/SciLifeLab/Sarek), developed at [SciLifeLab](https://github.com/SciLifeLab).

## Pipeline Flowchart
<p align="center">
  <img id="diagram" src="./docs/pipeline-flowchart.png"/>
</p>

## Directed Acyclic Graph
<p align="center">
  <img id="dag" src="./docs/dag.png"/>
</p>

##
<p align="center">
  <img src="./docs/brandenburg5_allegro.jpg">
</p>
