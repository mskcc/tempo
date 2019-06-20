FROM nfcore/base:latest

LABEL \
  authors="Yixiao Gong (gongy@mskcc.org)" 

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a

ENV PATH /opt/conda/envs/strelka2-manta-bcftools-vt/bin:$PATH

COPY filter-vcf.py /opt/conda/envs/strelka2-manta-bcftools-vt/bin
RUN chmod +x /opt/conda/envs/strelka2-manta-bcftools-vt/bin/filter-vcf.py
