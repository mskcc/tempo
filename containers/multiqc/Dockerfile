FROM ewels/multiqc:1.9

LABEL \
  authors="Anne Marie Noronha (noronhaa@mskcc.org)" \
  version.image="0.0.7" \
  version.general_stats_parse="0.0.2"

RUN apk update 
RUN echo "http://dl-8.alpinelinux.org/alpine/edge/community" >> /etc/apk/repositories 
RUN apk add gcc g++ gfortran 
RUN apk add python3-dev
RUN wget https://bootstrap.pypa.io/get-pip.py && python get-pip.py
RUN pip install --no-cache-dir pandas
RUN apk del --purge build-base libgfortran libpng-dev freetype-dev python3-dev py-pip

RUN mkdir /usr/bin/multiqc_custom_config
COPY wgs_multiqc_config.yaml /usr/bin/multiqc_custom_config
COPY exome_multiqc_config.yaml /usr/bin/multiqc_custom_config
COPY tempoLogo.png /usr/bin/multiqc_custom_config
RUN chmod 744 /usr/bin/multiqc_custom_config/*
COPY general_stats_parse.py /usr/bin/general_stats_parse.py
RUN chmod 755 /usr/bin/general_stats_parse.py
