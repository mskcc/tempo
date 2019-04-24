FROM broadinstitute/gatk:4.1.1.0

LABEL authors="Nikhil Kumar (kumarn1@mskcc.org), Evan Biederstedt (evan.biederstedt@gmail.com)" \
      version.image="1.0.0" \
      version.conpair="0.3.3" \
      version.gatk="4.1.1.0" \
      source.conpair="https://github.com/mskcc/Conpair/releases/tag/0.3.3" \
      source.r="https://pkgs.alpinelinux.org/package/edge/community/x86/R"      

ENV CONPAIR_VERSION 0.3.3

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
