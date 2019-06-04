FROM halllab/bcftools:v1.9

LABEL maintainer="C. Allan Bolipata (bolipatc@mskcc.org)" \
    version.image="1.0.0" \
    version.vt="0.57721" \
    version.filter-vcf="0.1.0"

RUN apt-get update && \
    apt-get install --yes \
    procps \
    gcc \
    make \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libssl-dev \
    libcurl4-openssl-dev \
    g++ \
    git \
    python-pip && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

RUN pip install pysam

RUN git clone https://github.com/atks/vt.git && \
    cd vt && \
    make && \
    cp /vt/vt /usr/bin

COPY filter-vcf.py /usr/bin
RUN chmod +x /usr/bin/filter-vcf.py