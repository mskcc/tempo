FROM continuumio/miniconda3:4.8.2

LABEL maintainer="Anne Marie Noronha (noronhaa@mskcc.org)" \
    version.base="4.8.2" \
    version.image="0.0.1" \
    version.samtools="1.9" \
    version.gatk="4.1.9.0"

ENV envName "gatk-samtools"

RUN apt-get update && apt-get install -y procps && apt-get clean -y

COPY environment.yml /
RUN conda update -n base -c defaults conda
RUN conda env create --name ${envName} -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/${envName}/bin:$PATH
RUN echo "export PATH=$PATH" > /etc/environment

