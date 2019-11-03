FROM python:3.7

LABEL \
    maintainer="Philip Jonsson (jonssonp@mskcc.org)" \
    contributor="C. Allan Bolipata (bolipatc@mskcc.org)" \
    version.image="1.0.0" \
    version.vt="0.57721" \
    version.filter-vcf="0.2.0" \
    version.fix-manta-vcf="0.1.0" \
    version.get-sv-info="0.1.0" \
    version.pysam="0.15.3" \
    source.getBaseCountsMultiSample="https://github.com/zengzheng123/GetBaseCountsMultiSample/releases/tag/v1.2.2" \
    version.getBaseCountsMultiSample="1.2.2" \
    version.survivor="1.0.7"

ENV GBCMS_VERSION 1.2.2
ENV VT_VERSION 0.57721

# Dependencies and tools
RUN apt-get update \
    && apt-get install --assume-yes cmake bedtools bcftools tabix

RUN pip install Cython pysam sklearn

# Install svtk
RUN git clone https://github.com/talkowski-lab/svtk.git \
    && cd svtk \
    && pip install -e .

# Install vt
RUN wget https://github.com/atks/vt/archive/${VT_VERSION}.tar.gz \
    && tar -xzf ${VT_VERSION}.tar.gz \
    && cd vt-${VT_VERSION} \
    && make \
    && cp vt /usr/bin

# Install GetBaseCountsMultiSample, copied from https://github.com/mskcc/roslin-variant/blob/2.5.x/build/containers/getbasecountsmultisample/1.2.2/Dockerfile
RUN cd /tmp \
    && wget https://github.com/zengzheng123/GetBaseCountsMultiSample/archive/v${GBCMS_VERSION}.zip \
    && unzip v${GBCMS_VERSION}.zip \
    # Install bamtools
    && cd /tmp/GetBaseCountsMultiSample-${GBCMS_VERSION}/bamtools-master \
    && rm -r build/ \
    && mkdir build \
    && cd build/ \
    && cmake -DCMAKE_CXX_FLAGS=-std=c++03 .. \
    && make \
    && make install \
    && cp ../lib/libbamtools.so.2.3.0 /usr/lib/ \
    # Install GetBaseCountsMultiSample itself
    && cd /tmp/GetBaseCountsMultiSample-${GBCMS_VERSION} \
    && make \
    && cp GetBaseCountsMultiSample /usr/bin/ 

# Add custom scripts
COPY filter-vcf.py /usr/bin
COPY fix-manta-vcf.py /usr/bin
COPY get-sv-info.py /usr/bin
RUN chmod +x /usr/bin/filter-vcf.py
RUN chmod +x /usr/bin/fix-manta-vcf.py
RUN chmod +x /usr/bin/get-sv-info.py
