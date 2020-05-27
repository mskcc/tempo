FROM ewels/multiqc:1.8

LABEL \
  authors="Anne Marie Noronha (noronhaa@mskcc.org)" \
  version.image="0.0.1"


RUN mkdir /usr/bin/multiqc_custom_config
COPY wgs_multiqc_config.yaml /usr/bin/multiqc_custom_config
COPY exome_multiqc_config.yaml /usr/bin/multiqc_custom_config
