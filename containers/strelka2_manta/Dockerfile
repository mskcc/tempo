FROM centos:6

MAINTAINER Jason Hwee <hweej@mskcc.org>

RUN yum -y update && \
    yum -y install wget && \
    yum -y install tar.x86_64 && \
    yum clean all

RUN wget https://github.com/Illumina/strelka/releases/download/v2.9.10/strelka-2.9.10.centos6_x86_64.tar.bz2 \
    && tar xvjf strelka-2.9.10.centos6_x86_64.tar.bz2 \
    && wget --no-check-certificate https://github.com/Illumina/manta/releases/download/v1.5.0/manta-1.5.0.centos6_x86_64.tar.bz2 \
    && tar xvjf manta-1.5.0.centos6_x86_64.tar.bz2

ENV PATH="/strelka-2.9.10.centos6_x86_64/bin:/manta-1.5.0.centos6_x86_64/bin:${PATH}"
ENV STRELKA_INSTALL_PATH="/strelka-2.9.10.centos6_x86_64"
