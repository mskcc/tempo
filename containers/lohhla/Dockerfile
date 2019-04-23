FROM continuumio/anaconda3:5.3.0

LABEL authors="Evan Biederstedt (evan.biederstedt@gmail.com)" \
      version.image="1.0.0" \
      version.lohhla"1.0.0"   

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
RUN conda install lohhla
RUN conda update lohhla


