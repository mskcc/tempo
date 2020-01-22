FROM cmopipeline/lohhla:1.1.3

LABEL authors="Yixiao Gong (gongy@mskcc.com)" \
      version.image="1.1.4" \
      version.lohhla="1.1.4" \
      source.lohhla="https://github.com/mskcc/lohhla/releases/tag/1.1.4"

WORKDIR /

RUN rm -rf /lohhla

RUN wget https://github.com/mskcc/lohhla/archive/1.1.4.zip && unzip 1.1.4.zip && mv /lohhla-1.1.4 /lohhla

WORKDIR /lohhla
