FROM halllab/svtools:v0.5.1

LABEL maintainer="Anne Marie Noronha (noronhaa@mskcc.org)" \
    version.image="0.0.3" 

RUN apt-get update && apt-get install -y procps gcc libz-dev
RUN apt-get install bedtools
RUN apt-get install -y sqlite3
RUN pip install --upgrade pip
RUN pip install pybedtools

ENV ANNOTSV /opt/annotsv

RUN mkdir -p /tmp \
    && apt-get update && \
    apt-get install --yes --allow-unauthenticated \
    bedtools \
    tcl \
    curl \
    wget \
    build-essential \
    bcftools && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

ENV ANNOTSV_VERSION v3.1.1

RUN mkdir -p unpack_annotsv \
    && wget https://github.com/lgmgeo/AnnotSV/archive/refs/tags/${ANNOTSV_VERSION}.tar.gz \
    && tar -C unpack_annotsv -xvzf ${ANNOTSV_VERSION}.tar.gz \
    && mv unpack_annotsv/AnnotSV* ${ANNOTSV} \
    && rm -rf unpack_annotsv \
    && cd ${ANNOTSV} \
    && make PREFIX=. install

