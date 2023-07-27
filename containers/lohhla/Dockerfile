FROM cmopipeline/lohhla:1.1.3

LABEL authors="Yixiao Gong (gongy@mskcc.com)" \
      version.image="1.1.7" \
      version.lohhla="1.1.7" \
      source.lohhla="https://github.com/mskcc/lohhla/releases/tag/1.1.7"

ENV LOHHLA_VERSION 1.1.7

RUN R -e "install.packages(c('plyr'), repos = 'http://cran.r-project.org')"

WORKDIR /

RUN rm -rf /lohhla

RUN wget https://github.com/mskcc/lohhla/archive/${LOHHLA_VERSION}.zip && unzip ${LOHHLA_VERSION}.zip && mv /lohhla-${LOHHLA_VERSION} /lohhla

WORKDIR /lohhla
