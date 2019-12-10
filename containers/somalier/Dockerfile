FROM ubuntu:14.04

LABEL maintainer="Yixiao Gong (gong@mskcc.org)" \
      version.ubuntu="14.04" \
      version.somalier="v0.2.6"

# install required packages
RUN apt-get update && apt-get install -y \
    build-essential \
    cmake \
    g++ \
    gfortran \
    zip \
    wget \
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

# install alfred
RUN cd /opt \
    && wget https://github.com/brentp/somalier/releases/download/v0.2.6/somalier \
    && chmod +x /opt/somalier

ENV PATH="/opt/:${PATH}"
