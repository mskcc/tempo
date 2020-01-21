FROM nfcore/base:latest

LABEL \
  authors="Yixiao Gong (gongy@mskcc.org)" \
  version.image="1.1.4" \
  version.lohhla="1.1.4" \
  source.lohhla="https://github.com/mskcc/lohhla/releases/tag/1.1.4"

ENV PATH /opt/conda/envs/env/bin:$PATH

## Set up channels
## https://bioconda.github.io/index.html
RUN conda config --add channels bioconda
RUN conda config --add channels conda-forge

## install unzip
RUN conda install -c conda-forge unzip

## we do this solely to make installing the LOHHLA dependencies easier
RUN conda install lohhla
RUN conda update lohhla
### the above installs all LOHHLA dependencies, e.g. samtools, picard, bedtools, Novalign, etc.

RUN conda clean -a

## download the mskcc repo of lohhla to run on Nextflow
## only takes single TN pair
## no working directory
## RUN git clone --single-branch --branch feature/work_with_Nexflow https://github.com/mskcc/lohhla.git

RUN wget https://github.com/mskcc/lohhla/archive/1.1.4.zip && unzip 1.1.4.zip && mv /lohhla-1.1.4 /lohhla

RUN pwd

WORKDIR /lohhla
