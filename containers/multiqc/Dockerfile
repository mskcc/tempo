FROM ewels/multiqc:1.9

LABEL \
  authors="Anne Marie Noronha (noronhaa@mskcc.org)" \
  version.image="0.0.4"

RUN mkdir /usr/bin/multiqc_custom_config
COPY wgs_multiqc_config.yaml /usr/bin/multiqc_custom_config
COPY exome_multiqc_config.yaml /usr/bin/multiqc_custom_config
COPY tempoLogo.png /usr/bin/multiqc_custom_config
RUN chmod 544 /usr/bin/multiqc_custom_config/*
