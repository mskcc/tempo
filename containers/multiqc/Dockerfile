FROM continuumio/miniconda3:4.10.3p0

LABEL maintainer="Anne Marie Noronha (noronhaa@mskcc.org)" \
    version.base="4.8.2" \
    version.image="0.1.3" \
    version.general_stats_parse="0.0.4" \
    version.parse_alfred="0.0.2" \
    version.clean_fastp="0.0.1"

ENV envName "multiqc"
ENV MULTIQC_VERSION "1.9"

RUN apt-get update && apt-get install -y procps && apt-get clean -y

COPY environment.yml /
RUN conda update -n base -c defaults conda
RUN conda env create --name ${envName} -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/${envName}/bin:$PATH
RUN echo "export PATH=$PATH" > /etc/environment

RUN mkdir -p /tmp/ \
  && cd /tmp \
  && wget https://github.com/ewels/MultiQC/archive/v${MULTIQC_VERSION}.tar.gz \
  && tar -xzf v${MULTIQC_VERSION}.tar.gz \
  && mkdir MultiQC-${MULTIQC_VERSION}/multiqc/templates/tempo \
  && mkdir -p MultiQC-${MULTIQC_VERSION}/multiqc/templates/tempo/assets/css/ \
  && cp MultiQC-${MULTIQC_VERSION}/multiqc/templates/default/__init__.py MultiQC-${MULTIQC_VERSION}/multiqc/templates/tempo/ \
  && echo "template_parent = 'default'" >> MultiQC-${MULTIQC_VERSION}/multiqc/templates/tempo/__init__.py  
COPY default_multiqc.css /tmp/MultiQC-${MULTIQC_VERSION}/multiqc/templates/tempo/assets/css/
COPY setup.py /tmp/MultiQC-${MULTIQC_VERSION}/
RUN cd /tmp/MultiQC-${MULTIQC_VERSION}/ && pip install .
ENV configPath /opt/conda/envs/multiqc/lib/python3.8/site-packages/multiqc/example_tempo_configs
RUN mkdir -p ${configPath}
COPY example_wgs_multiqc_config.yaml ${configPath}
COPY example_exome_multiqc_config.yaml ${configPath}
RUN chmod 744 ${configPath}/*

# copy auxiliary scripts
COPY general_stats_parse.py /opt/conda/envs/${envName}/bin/general_stats_parse.py
COPY parse_alfred.py /opt/conda/envs/${envName}/bin/parse_alfred.py
COPY clean_fastp.py /opt/conda/envs/${envName}/bin/clean_fastp.py
RUN cd /opt/conda/envs/${envName}/bin/ && chmod 755 general_stats_parse.py parse_alfred.py clean_fastp.py

