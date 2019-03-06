FROM alpine:3.8

LABEL maintainer="Nikhil Kumar (kumarn1@mskcc.org)" \
      version.image="1.0.0" \
      version.delly="0.7.7" \
      version.htslib="1.6" \
      version.bcftools="1.6" \
      version.samtools="1.6" \
      version.alpine="3.8" \
      source.delly="https://github.com/dellytools/delly/releases/tag/v0.7.7" \
      source.htslib="https://github.com/samtools/htslib/releases/tag/1.6" \
      source.bcftools="https://github.com/samtools/bcftools/releases/tag/1.6" \
      source.samtools="https://github.com/samtools/samtools/releases/tag/1.6"

ENV DELLY_VERSION v0.8.1


# set environment
ENV BOOST_ROOT /usr
ENV SEQTK_ROOT /opt/htslib

RUN apk add --update --no-cache ncurses \
    && apk add --virtual=deps --update --no-cache ncurses-dev musl-dev git g++ make zlib-dev \
    && apk add bash build-base asciidoc bzip2-dev xz-dev boost-dev \
    && mkdir /opt \
    && cd /opt \
    # install delly
    && cd /opt \
    && git clone --recursive https://github.com/dellytools/delly.git -b ${DELLY_VERSION} \
    && cd /opt/delly/ \
    && make all \
    && install -p /opt/delly/src/delly /usr/local/bin/ 

