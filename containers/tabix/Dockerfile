FROM ubuntu:18.04

LABEL maintainer="C. Allan Bolipata (bolipatc@mskcc.org)" \
      version.image="1.0.0" \
      version.htstools="0.1.1" \
      version.htslib="1.5" \
      version.ubuntu="18.04" \
      source.htstools="https://github.com/mskcc/htstools/releases/tag/snp_pileup_0.1.1" \
      source.htslib="https://github.com/samtools/htslib/releases/tag/1.5"

RUN apt-get update \
    && apt-get -y upgrade \
      # install all the build-related tools
            && apt-get install -y build-essential \
            && apt-get install -y software-properties-common \
            && apt-get install -y byobu curl git htop man unzip wget \
      # install tabix
            && apt-get install -y tabix


