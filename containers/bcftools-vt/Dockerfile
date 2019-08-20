FROM halllab/bcftools:v1.9

LABEL maintainer="C. Allan Bolipata (bolipatc@mskcc.org)" \
    contributor="Philip Jonsson (jonssonp@mskcc.org)" \
    version.image="1.2.0" \
    version.vt="0.57721" \
    version.filter-vcf="0.2.0" \
    version.pysam="0.15.2" \
    source.getBaseCountsMultiSample="https://github.com/zengzheng123/GetBaseCountsMultiSample/releases/tag/v1.2.2" \
    version.getBaseCountsMultiSample="1.2.2"

ENV GBCMS_VERSION 1.2.2

RUN apt-get update && \
    apt-get install --yes \
    procps \
    gcc \
    make \
    cmake \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libssl-dev \
    libcurl4-openssl-dev \
    g++ \
    git \
    wget \
    zip \
    python-pip && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Install pysam
RUN pip install pysam==0.15.2

# Install vt
RUN git clone https://github.com/atks/vt.git && \
    cd vt && \
    make && \
    cp /vt/vt /usr/bin

# Install GetBaseCountsMultiSample, copied from https://github.com/mskcc/roslin-variant/blob/2.5.x/build/containers/getbasecountsmultisample/1.2.2/Dockerfile
RUN cd /tmp && \
    wget https://github.com/zengzheng123/GetBaseCountsMultiSample/archive/v${GBCMS_VERSION}.zip && \
    unzip v${GBCMS_VERSION}.zip && \
    # Install bamtools
    cd /tmp/GetBaseCountsMultiSample-${GBCMS_VERSION}/bamtools-master && \
    rm -r build/ && \
    mkdir build && \
    cd build/ && \
    cmake -DCMAKE_CXX_FLAGS=-std=c++03 .. && \
    make && \
    make install && \
    cp ../lib/libbamtools.so.2.3.0 /usr/lib/ && \
    # Install GetBaseCountsMultiSample itself
    cd /tmp/GetBaseCountsMultiSample-${GBCMS_VERSION} && \
    make && \
    cp GetBaseCountsMultiSample /usr/bin/ 

# Add filter script
COPY filter-vcf.py /usr/bin
RUN chmod +x /usr/bin/filter-vcf.py
