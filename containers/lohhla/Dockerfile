FROM continuumio/anaconda3:5.3.0

LABEL authors="Evan Biederstedt (evan.biederstedt@gmail.com)" \
      version.image="1.0.0" \
      version.lohhla="1.0.0" \
      source.msig="https://github.com/mskcc/lohhla/releases/tag/1.0.0"

RUN conda create -n env python=3.6
RUN echo "source activate env" > ~/.bashrc

ENV PATH /opt/conda/envs/env/bin:$PATH

# Updating Anaconda packages
RUN conda update conda 
RUN conda update anaconda
RUN conda update --all

## Set up channels
## https://bioconda.github.io/index.html
RUN conda config --add channels defaults
RUN conda config --add channels bioconda
RUN conda config --add channels conda-forge

## install lohhla
## https://bioconda.github.io/recipes/lohhla/README.html

## we do this solely to make installing the LOHHLA dependencies easier

RUN conda install lohhla
RUN conda update lohhla


### the above installs all LOHHLA dependencies, e.g. samtools, picard, bedtools, Novalign, etc. 

RUN rm /lohhla

## download the mskcc repo of lohhla to run on Nextflow
## only takes single TN pair
## no working directory
## RUN git clone --single-branch --branch feature/work_with_Nexflow https://github.com/mskcc/lohhla.git
RUN wget https://github.com/mskcc/lohhla/archive/1.0.0.zip && mv lohhla-1.0.0 lohhla 
