FROM continuumio/miniconda3:4.8.2

LABEL maintainer="Anne Marie Noronha (noronhaa@mskcc.org)" \
    version.base="4.8.2" \
    version.image="0.1.2" \
    version.general_stats_parse="0.0.4" \
    version.parse_alfred="0.0.2"

ENV envName "multiqc"

RUN apt-get update && apt-get install -y procps && apt-get clean -y

COPY environment.yml /
RUN conda update -n base -c defaults conda
RUN conda env create --name ${envName} -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/${envName}/bin:$PATH
RUN echo "export PATH=$PATH" > /etc/environment

RUN mkdir -p /tmp/ \
  && cd /tmp \
  && wget https://github.com/ewels/MultiQC/archive/v1.9.tar.gz \
  && tar -xzf v1.9.tar.gz \
  && mkdir MultiQC-1.9/multiqc/templates/tempo \
  && mkdir -p MultiQC-1.9/multiqc/templates/tempo/assets/css/ \
  && cp MultiQC-1.9/multiqc/templates/default/__init__.py MultiQC-1.9/multiqc/templates/tempo/ \
  && echo "template_parent = 'default'" >> MultiQC-1.9/multiqc/templates/tempo/__init__.py  
COPY default_multiqc.css /tmp/MultiQC-1.9/multiqc/templates/tempo/assets/css/
COPY setup.py /tmp/MultiQC-1.9/
RUN cd /tmp/MultiQC-1.9/ && pip install .
ENV configPath /opt/conda/envs/multiqc/lib/python3.8/site-packages/multiqc/example_tempo_configs
RUN mkdir -p ${configPath}
COPY example_wgs_multiqc_config.yaml ${configPath}
COPY example_exome_multiqc_config.yaml ${configPath}
RUN chmod 744 ${configPath}/*
COPY general_stats_parse.py /opt/conda/envs/${envName}/bin/general_stats_parse.py
RUN chmod 755 /opt/conda/envs/${envName}/bin/general_stats_parse.py
COPY parse_alfred.py /opt/conda/envs/${envName}/bin/parse_alfred.py
RUN chmod 755 /opt/conda/envs/${envName}/bin/parse_alfred.py

