FROM nfcore/base:latest

LABEL \
  authors="Yixiao Gong (gongy@mskcc.org)" \
  contributors="Anne Marie Noronha (noronhaa@mskcc.org)" \
  version.image="2.0.1"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a

ENV PATH /opt/conda/envs/strelka2-manta-bcftools-vt/bin:$PATH
