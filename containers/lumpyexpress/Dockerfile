FROM ubuntu:18.04
LABEL maintainer "C. Allan Bolipata <bolipatc@mskcc.org>"

ENV LUMPYSV_VERSION 0.3.0
ENV SAMBLASTER_VERSION v.0.1.24
ENV SAMTOOLS_VERSION 1.9
ENV SAMBAMBA_VERSION 0.6.8

# Build dependencies
RUN apt-get update -qq \
    && apt-get -y install \
        zlib1g-dev bzip2 libbz2-dev liblzma-dev libncurses5-dev \
        python python-pip python-dev python3 wget ldc \
        gawk build-essential \
        libcurl4-gnutls-dev \
        autoconf \
        libssl-dev \
        git \
    && pip install --upgrade pip numpy scipy pysam

# install samtools
RUN cd /tmp \
    && wget https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 \
    && tar vxjf samtools-${SAMTOOLS_VERSION}.tar.bz2 \
    && cd samtools-${SAMTOOLS_VERSION} \
    && ./configure --prefix=/usr/bin \
    && make \
    && rm -rf /tmp/*

# Install sambamba
 RUN cd /tmp \
    && wget https://github.com/biod/sambamba/releases/download/v0.6.8/sambamba-0.6.8-linux-static.gz \
    && gunzip sambamba-0.6.8-linux-static.gz \
    && mv sambamba-0.6.8-linux-static /usr/bin/sambamba \
    && chmod +x /usr/bin/sambamba

# Install samblaster
 RUN cd /tmp \
    && wget https://github.com/GregoryFaust/samblaster/archive/${SAMBLASTER_VERSION}.tar.gz \
    && tar xvzf ${SAMBLASTER_VERSION}.tar.gz \
    && cd /tmp/samblaster-${SAMBLASTER_VERSION} \
    && make \
    && cp samblaster /usr/bin \
    && rm -rf /tmp/*

# Install LUMPY
RUN cd /opt \
    && git clone --recursive https://github.com/arq5x/lumpy-sv.git \
    && cd lumpy-sv/ \
    && git checkout tags/${LUMPYSV_VERSION} \
    && make \
    && cd bin/ \
    # set lumpyexpress.config
    && sed -i 's/SAMTOOLS=/SAMTOOLS=\/usr\/bin\/samtools/' lumpyexpress.config \
    && sed -i 's/SAMBLASTER=/SAMBLASTER=\/usr\/bin\/samblaster/' lumpyexpress.config \
    && sed -i 's/SAMBAMBA=/SAMBAMBA=\/usr\/bin\/sambamba/' lumpyexpress.config

env PATH=/opt/lumpy-sv/bin:$PATH
