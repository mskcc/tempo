FROM continuumio/miniconda3:4.8.2

LABEL maintainer="Anne Marie Noronha (noronhaa@mskcc.org)" \
    version.base="4.8.2" \
    version.image="0.0.9" \
    version.general_stats_parse="0.0.3"

ENV envName "multiqc"

RUN apt-get update && apt-get install -y procps && apt-get clean -y

COPY environment.yml /
RUN conda update -n base -c defaults conda
RUN conda env create --name ${envName} -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/${envName}/bin:$PATH

ENV configPath /opt/conda/envs/multiqc/lib/python3.8/site-packages/multiqc/example_tempo_configs
RUN mkdir -p ${configPath}
COPY example_wgs_multiqc_config.yaml ${configPath}
COPY example_exome_multiqc_config.yaml ${configPath}
RUN chmod 744 ${configPath}/*
COPY general_stats_parse.py /opt/conda/envs/${envName}/bin/general_stats_parse.py
RUN chmod 755 /opt/conda/envs/${envName}/bin/general_stats_parse.py
