FROM ubuntu:18.04

LABEL maintainer="C. Allan Bolipata (bolipatc@mskcc.org)" \
      version.ubuntu="18.04" \
      version.tabix="1.7-2"

RUN apt-get update \
    && apt-get -y upgrade \
      # install all the build-related tools
            && apt-get install -y build-essential \
            && apt-get install -y software-properties-common \
            && apt-get install -y byobu curl git htop man unzip wget \
      # install tabix
            && apt-get install -y tabix


