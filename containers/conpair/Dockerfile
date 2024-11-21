FROM ubuntu:16.04

LABEL authors="Nikhil Kumar (kumarn1@mskcc.org), Evan Biederstedt (evan.biederstedt@gmail.com), C. Allan Bolipata (bolipatc@mskcc.org)" \
      version.image="1.0.0" \
      version.conpair="0.3.3" \
      version.gatk="3.8-1" \
      source.conpair="https://github.com/mskcc/Conpair/releases/tag/0.3.3" \
      source.r="https://pkgs.alpinelinux.org/package/edge/community/x86/R"      

ENV CONPAIR_VERSION 0.3.3

RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get install -y python && \
    apt-get install -y wget \
    curl \
    bc \
    unzip \
    less \
    bedtools \
    samtools \
    openjdk-8-jdk \
    tabix \
    bzip2 \
    software-properties-common && \
    apt-get -y clean  && \
    apt-get -y autoclean  && \
    apt-get -y autoremove

# Download GATK
RUN cd /tmp \ 
    && wget "https://software.broadinstitute.org/gatk/download/auth?package=GATK-archive&version=3.8-1-0-gf15c1c3ef" -O gatk-3.8-1.tar.bz2 \ 
    && tar xvjf gatk-3.8-1.tar.bz2 \
    && mv /tmp/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar /usr/bin/GenomeAnalysisTK.jar \
    && rm -rf /tmp/*

## install R 
RUN apt-get install -y r-base-core r-base-dev \
    ## install R dependencies
    && R -e "install.packages('ggplot2', dependencies=TRUE, repos='http://cran.rstudio.com/')" \
    && R -e "install.packages('reshape2', dependencies=TRUE, repos='http://cran.rstudio.com/')" 

## Python dependencies and Java 
## install numpy and scip
RUN pip install numpy==1.15.4 \
    && pip install scipy==1.1.0

## download Conpair
RUN cd /tmp && wget https://github.com/mskcc/Conpair/archive/${CONPAIR_VERSION}.tar.gz \
    && tar xvzf ${CONPAIR_VERSION}.tar.gz \
    ## install conpair
    && mv /tmp/Conpair-${CONPAIR_VERSION} /usr/bin/conpair \
    ## clean up
    && rm -rf /tmp/*

ENV PYTHONNOUSERSITE set
