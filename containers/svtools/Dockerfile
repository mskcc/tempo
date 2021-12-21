FROM halllab/svtools:v0.5.1

LABEL maintainer="Anne Marie Noronha (noronhaa@mskcc.org)" \
    version.image="0.0.3" 

RUN mkdir -p /tmp \
    && apt-get update \ 
    && apt-get install -y \
    procps \
    gcc \
    libz-dev \
    sqlite3 \
    bedtools \
    tcl \
    curl \
    wget \
    build-essential \
    bcftools \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

RUN pip install --upgrade pip
RUN pip install pybedtools

ENV ANNOTSV /opt/annotsv
ENV ANNOTSV_VERSION v3.1.1

# This is an installation of AnnotSV source code, with the preparation of the reference files intentionally omitted
# Reference files were prepared outside of docker using this container after it was built.
# To prepare reference files, start a container and: cd /opt/annotsv/ && make PREFIX=<install> install-human-annotation
RUN mkdir -p unpack_annotsv \
    && wget https://github.com/lgmgeo/AnnotSV/archive/refs/tags/${ANNOTSV_VERSION}.tar.gz \
    && tar -C unpack_annotsv -xvzf ${ANNOTSV_VERSION}.tar.gz \
    && mv unpack_annotsv/AnnotSV* ${ANNOTSV} \
    && rm -rf unpack_annotsv \
    && cd ${ANNOTSV} \
    && make PREFIX=. install

