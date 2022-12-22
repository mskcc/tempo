FROM continuumio/miniconda3:4.12.0

LABEL maintainer="Anne Marie Noronha (noronhaa@mskcc.org)" \
    version.base="4.10.3" \
    version.image="0.0.1" 

ENV envName "svclone"

RUN apt-get update && apt-get install -y procps && apt-get clean -y

RUN conda update -n base -c defaults conda
RUN conda create --name ${envName} -c bioconda -c conda-forge svclone=1.1.1-0 r-base
ENV PATH /opt/conda/envs/${envName}/bin:$PATH
RUN echo "export PATH=$PATH" > /etc/environment

RUN mkdir -p /config && \
    cd /config && \
    wget https://github.com/mcmero/SVclone/raw/4a02d5e2d8f548b2f38e9f7f77ea432bf456b349/svclone_config.ini 

RUN sed -i "s/itrx_class: INTRX/itrx_class: INTRX,TRA,BND/g" /config/svclone_config.ini

# The following solves a bug where certain files are written to the wrong location. Once svclone is upgraded from 1.1.1-0, this line should be deleted.
RUN wget https://raw.githubusercontent.com/mcmero/SVclone/034d04996216f4a3378f23e583c984724e05ffbc/SVclone/post_assign.R -O /opt/conda/envs/${envName}/lib/python3.6/site-packages/SVclone/post_assign.R
