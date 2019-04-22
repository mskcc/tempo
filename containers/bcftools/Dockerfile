FROM halllab/bcftools:v1.9

MAINTAINER C. Allan Bolipata <bolipatc@mskcc.org>

RUN apt-get update && \
    apt-get install --yes \
    procps && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*
