FROM cmopipeline/lohhla:1.1.3

LABEL authors="Yixiao Gong (gongy@mskcc.com)" \
      version.image="1.1.6" \
      version.lohhla="1.1.6" \
      source.lohhla="https://github.com/mskcc/lohhla/releases/tag/1.1.6"

RUN R -e "install.packages(c('plyr'), repos = 'http://cran.r-project.org')"

WORKDIR /

RUN rm -rf /lohhla

RUN wget https://github.com/mskcc/lohhla/archive/1.1.6.zip && unzip 1.1.6.zip && mv /lohhla-1.1.6 /lohhla

WORKDIR /lohhla
