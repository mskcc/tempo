FROM nfcore/base:latest

LABEL \
  authors="Yixiao Gong (gongy@mskcc.org)" 

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a

ENV PATH /opt/conda/envs/fastp-bwa-samtools/bin:$PATH
