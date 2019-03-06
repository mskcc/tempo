FROM ubuntu:14.04

LABEL maintainer="C. Allan Bolipata (bolipatc@mskcc.org)" \
      version.image="1.0.0" \
      version.htstools="0.1.1" \
      version.htslib="1.5" \
      version.ubuntu="14.04" \
      source.htstools="https://github.com/mskcc/htstools/releases/tag/snp_pileup_0.1.1" \
      source.htslib="https://github.com/samtools/htslib/releases/tag/1.5"

ENV HTSTOOLS_VERSION 0.1.1
ENV HTSLIB_VERSION 1.5


RUN apt-get update \
    && apt-get -y upgrade \
      # install all the build-related tools
            && apt-get install -y build-essential \
            && apt-get install -y software-properties-common \
            && apt-get install -y byobu curl git htop man unzip wget \
      # install system packages that the Perl modules will need
            && apt-get install -y zlib1g-dev bzip2 libbz2-dev liblzma-dev \
      # download, unzip, install htslib
            && cd /tmp \
            && wget https://github.com/samtools/htslib/releases/download/${HTSLIB_VERSION}/htslib-${HTSLIB_VERSION}.tar.bz2 \
            && tar xvjf htslib-${HTSLIB_VERSION}.tar.bz2 \
            && cd /tmp/htslib-${HTSLIB_VERSION} \
            && ./configure \
            && make && make install \
      # download, unzip
            && cd /tmp \
            && wget https://github.com/mskcc/htstools/archive/snp_pileup_${HTSTOOLS_VERSION}.tar.gz \
            && tar xvzf snp_pileup_${HTSTOOLS_VERSION}.tar.gz \
            && cd /tmp/htstools-snp_pileup_${HTSTOOLS_VERSION} \
            && cp /tmp/htslib-${HTSLIB_VERSION}/libhts.so /usr/lib \
            && cp /tmp/htslib-${HTSLIB_VERSION}/libhts.so.2 /usr/lib \
      # install snp-pileup and ppflag-fixer
            && cd /tmp/htstools-snp_pileup_${HTSTOOLS_VERSION} \
            && g++ -std=c++11 snp-pileup.cpp -lhts -o snp-pileup \
            && g++ -std=c++11 ppflag-fixer.cpp -lhts -o ppflag-fixer \
      #move executables into bin
            && cp /tmp/htstools-snp_pileup_${HTSTOOLS_VERSION}/snp-pileup /usr/bin \
            && cp /tmp/htstools-snp_pileup_${HTSTOOLS_VERSION}/ppflag-fixer /usr/bin \
      # add clean up
            && rm -rf /tmp/*

ENV PYTHONNOUSERSITE set
