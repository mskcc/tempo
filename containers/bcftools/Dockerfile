FROM halllab/bcftools:v1.9

MAINTAINER C. Allan Bolipata <bolipatc@mskcc.org>

RUN apt-get update && \
    apt-get install --yes \
    libncurses5-dev vcftools libbz2-dev liblzma-dev procps && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*
