FROM rocker/r-ubuntu:18.04

LABEL maintainer="Jun Woo <wooh@mskcc.org>" \
      version.tempoSig="0.2.1" \
      version.ubuntu="18.04" \
      source.msig="https://github.com/mskcc/tempoSig.git"

RUN apt-get update && apt-get install -y \
    apt-utils \
    build-essential \
    git \
    libgsl23 \
    libgsl-dev \
    libcurl4-openssl-dev \
    libssl-dev 

RUN R -e "install.packages(c('gtools','argparse','coneproj','data.table'), dependencies=TRUE)"

RUN git clone --branch docker --single-branch https://github.com/mskcc/tempoSig

RUN R -e "install.packages('tempoSig', repos=NULL, type='source')"

ENV PATH="/tempoSig/exec/:${PATH}"

WORKDIR /data
