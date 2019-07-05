FROM: r-base:3.5.3

RUN apt-get -y update && \
    apt-get -y install curl \ 
    openssl \
    libcurl4-openssl-dev \
    libssl-dev \
    procps

# Add BiocManager and install maftools
RUN R -e "source('https://bioconductor.org/biocLite.R'); biocLite('maftools')"


