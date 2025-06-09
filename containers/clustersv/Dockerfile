FROM continuumio/miniconda3:4.10.3

LABEL author="Anne Marie Noronha (noronhaa@mskcc.org)" \
    maintainer="Yixiao Gong (gongy@mskcc.org)" \
    version.base="4.8.2" \
    version.image="0.0.1"

RUN apt-get update -y && apt-get install -y procps

ENV envName "clusterSV"

COPY environment.yml /
RUN conda update -n base -c defaults conda
RUN conda env create --name ${envName} -f /environment.yml && conda clean -a
RUN apt-get install -y git
RUN cd /opt \
    && git clone https://github.com/cancerit/ClusterSV.git \
    && cd ClusterSV \
    && git reset --hard 1d7eeea8f22133d811bed6a4cbbaf1ad2122e2c4
ENV PATH /opt/conda/envs/${envName}/bin:$PATH

