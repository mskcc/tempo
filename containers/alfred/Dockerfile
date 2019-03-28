FROM ubuntu:14.04

LABEL maintainer="Sinisa Ivkovic (ivkovics@mskcc.org)" \
      version.ubuntu="14.04" \
      version.tabix="v0.1.17"

# install required packages
RUN apt-get update && apt-get install -y \
    build-essential \
    cmake \
    g++ \
    gfortran \
    git \
    hdf5-tools \
    libboost-date-time-dev \
    libboost-program-options-dev \
    libboost-system-dev \
    libboost-filesystem-dev \
    libboost-iostreams-dev \
    libbz2-dev \
    libhdf5-dev \
    libncurses-dev \
    liblzma-dev \
    zlib1g-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# set environment
ENV BOOST_ROOT /usr
ENV EBROOTHTSLIB /opt/htslib

# install alfred
RUN cd /opt \
    && git clone https://github.com/samtools/htslib.git \
    && cd /opt/htslib \
    && make \
    && make lib-static \
    && make install
RUN cd /opt \
    && git clone --branch v0.1.17 https://github.com/tobiasrausch/alfred.git \
    && cd /opt/alfred/ \
    && make STATIC=1 all \
    && make install

# install R packages
RUN sudo sh -c 'echo "deb http://cran.rstudio.com/bin/linux/ubuntu trusty/" >> /etc/apt/sources.list'
RUN gpg --keyserver keyserver.ubuntu.com --recv-key E084DAB9
RUN gpg -a --export E084DAB9 | sudo apt-key add -
RUN apt-get update && apt-get -y install r-base
RUN R -e "install.packages(c('ggplot2'), repos='http://cran.us.r-project.org')"

ENV PATH="/opt/alfred/bin:${PATH}"

