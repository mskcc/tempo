FROM alpine:3.8

LABEL authors="Nikhil Kumar (kumarn1@mskcc.org), Evan Biederstedt (evan.biederstedt@gmail.com)" \
      version.image="1.0.0" \
      version.conpair="v0.2" \
      version.alpine="3.8" \
      version.gatk="4.1.1.0" \
      source.gatk="https://github.com/broadinstitute/gatk/releases/download/4.1.1.0/gatk-4.1.1.0.zip" \
      source.conpair="https://github.com/nygenome/Conpair/releases/tag/v0.2" \
      source.r="https://pkgs.alpinelinux.org/package/edge/community/x86/R"      

ENV CONPAIR_VERSION 0.2
ENV GATK_VERSION 4.1.1.0


RUN apk add --update \
    # for wget
    && apk add ca-certificates openssl \
    # for compilation
    && apk add build-base musl-dev python py-pip python-dev py-setuptools\
    && cd /tmp 

## install R 
RUN apk add R R-dev \
    ## install R dependencies
    && R -e "install.packages('ggplot2', dependencies=TRUE, repos='http://cran.rstudio.com/')" \
    && R -e "install.packages('reshape2', dependencies=TRUE, repos='http://cran.rstudio.com/')" 

## Python dependencies and Java 
## install numpy and scipy
RUN pip install numpy==1.15.4 \
    && pip install scipy==1.1.0 \
    ## install java
    && apk add openjdk8-jre-base 

## download GATK 
RUN wget https://github.com/broadinstitute/gatk/releases/download/4.1.1.0/gatk-4.1.1.0.zip \
    && unzip gatk-4.1.1.0.zip \ 
    && cp gatk-${GATK_VERSION}/gatk-package-${GATK_VERSION}-local.jar /usr/bin/gatk.jar  \
    ## download and unzip conpair
    && cd /tmp && wget https://github.com/nygenome/Conpair/archive/v${CONPAIR_VERSION}.tar.gz \
    && tar xvzf ${CONPAIR_VERSION}.tar.gz \
    ## install conpair
    && mv /tmp/Conpair-${CONPAIR_VERSION} /usr/bin/conpair \
    ## clean up
    && rm -rf /tmp/*

ENV PYTHONNOUSERSITE set

ENTRYPOINT ["/bin/sh"]