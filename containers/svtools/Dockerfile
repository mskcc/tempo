FROM halllab/svtools:v0.5.1

LABEL maintainer="Anne Marie Noronha (noronhaa@mskcc.org)" \
    version.image="0.0.4" 

RUN mkdir -p /tmp \
    && apt-get update \
    && apt-get install -y \
    procps \
    gcc \
    libz-dev \
    sqlite3 \
    bedtools \
    curl \
    wget \
    build-essential \
    bcftools \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*


