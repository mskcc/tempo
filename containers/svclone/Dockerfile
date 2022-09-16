FROM continuumio/miniconda3:4.10.3

LABEL maintainer="Anne Marie Noronha (noronhaa@mskcc.org)" \
    version.base="4.10.3" \
    version.image="0.0.1" 

ENV envName "svclone"

RUN apt-get update && apt-get install -y procps && apt-get clean -y

COPY environment.yml /
RUN conda update -n base -c defaults conda
RUN conda env create --name ${envName} -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/${envName}/bin:$PATH
RUN echo "export PATH=$PATH" > /etc/environment

RUN mkdir -p /config && \
    cd /config && \
    wget https://github.com/mcmero/SVclone/raw/4a02d5e2d8f548b2f38e9f7f77ea432bf456b349/svclone_config.ini 

RUN sed -i "s/itrx_class: INTRX/itrx_class: INTRX,TRA,BND/g" /config/svclone_config.ini
