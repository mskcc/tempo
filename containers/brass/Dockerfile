FROM quay.io/wtsicgp/brass:v6.3.4

LABEL maintainer="Anne Marie Noronha (noronhaa@mskcc.org)" \
    version.image="0.0.1" 

USER root

RUN mkdir -p /tmp \
    && cd /tmp \
    && apt-get update \
    && apt-get install -y \
         procps \
	 bcftools \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

