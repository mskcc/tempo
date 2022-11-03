FROM quay.io/wtsicgp/brass:v6.3.4

LABEL maintainer="Anne Marie Noronha (noronhaa@mskcc.org)" \
    version.image="0.0.2" 

USER root

RUN mkdir -p /tmp \
    && cd /tmp \
    && apt-get update \
    && apt-get install -y \
         procps \
	 libz-dev \
         libbz2-dev liblzma-dev \
	 make \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

ENV HTSLIB_VERSION 1.9
ENV BCFTOOLS_VERSION 1.9

RUN cd /tmp && curl -L -o htslib-${HTSLIB_VERSION}.tar.bz2 https://github.com/samtools/htslib/releases/download/${HTSLIB_VERSION}/htslib-${HTSLIB_VERSION}.tar.bz2 \
            && tar xvjf htslib-${HTSLIB_VERSION}.tar.bz2 \
            && cd /tmp/htslib-${HTSLIB_VERSION} \
            && ./configure \
            && make && make install \
            && rm -rf /tmp/*

RUN cd /tmp && curl -L -o tmp2.tar.gz https://github.com/samtools/bcftools/releases/download/${BCFTOOLS_VERSION}/bcftools-${BCFTOOLS_VERSION}.tar.bz2 && \
    mkdir bcftools && \
    tar -C bcftools --strip-components 1 -jxf tmp2.tar.gz && \
    cd bcftools && \
    make && \
    make install && \
    cd .. && \
    rm -rf /tmp/*

