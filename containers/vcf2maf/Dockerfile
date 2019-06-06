FROM ensemblorg/ensembl-vep:release_96.0

MAINTAINER C. Allan Bolipata <bolipatc@mskcc.org>

ENV PATH $VEP_PATH/htslib:$PATH
ENV PERL5LIB $VEP_PATH:/opt/lib/perl5:$PERL5LIB

ENV SAMTOOLS_VERSION 1.9
ENV BCFTOOLS_VERSION 1.9
ENV VCF2MAF_VERSION 1.6.17
ENV VEP_VERSION 96.0

WORKDIR /tmp

USER root
RUN apt-get update

# Install samtools
RUN apt-get update && \
    apt-get install --yes \
    libncurses5-dev \
    vcftools \
    libbz2-dev \
    liblzma-dev && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

RUN curl -L -o tmp.tar.gz https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 && \
    mkdir samtools && \
    tar -C samtools --strip-components 1 -jxf tmp.tar.gz && \
    cd samtools && \
    ./configure && \
    make && \
    make install && \
    cd .. && \
    rm -rf *

RUN curl -L -o tmp2.tar.gz https://github.com/samtools/bcftools/releases/download/${BCFTOOLS_VERSION}/bcftools-${BCFTOOLS_VERSION}.tar.bz2 && \
    mkdir bcftools && \
    tar -C bcftools --strip-components 1 -jxf tmp2.tar.gz && \
    cd bcftools && \
    make && \
    make install && \
    cd .. && \
    rm -rf *

# Install vcf2maf
WORKDIR /opt/

RUN curl -ksSL -o tmp.tar.gz https://github.com/mskcc/vcf2maf/archive/v${VCF2MAF_VERSION}.tar.gz && \
    tar --strip-components 1 -zxf tmp.tar.gz && \
    rm tmp.tar.gz && \
    chmod +x *.pl

# Install R with depencies
RUN apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get install --yes \
    r-base-core r-base-dev && \
    R -e "install.packages('data.table', dependencies = TRUE, repos = 'http://cran.rstudio.com')"

# Add filter script
COPY filter-maf.R /usr/bin
RUN chmod +x /usr/bin/filter-maf.R
