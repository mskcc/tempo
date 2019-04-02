FROM centos:latest

LABEL maintainer="Jason Hwee <hweej@mskcc.org>" \
    version.image="1.0.0" \
    version.centos="7-latest" \
    version.fastp="0.19.7" \
    source.fastp="https://github.com/OpenGene/fastp.git"

ENV FASTP_VERSION_TAG v0.19.7

# Install Dependencies
RUN yum -y update \ 
    && yum install -y git gcc-c++ make zlib-devel \
    && yum clean all

# Fetch fastp code
RUN git clone --branch ${FASTP_VERSION_TAG} https://github.com/OpenGene/fastp.git

# Compile && Install
WORKDIR /fastp
RUN make && make install

WORKDIR /data
