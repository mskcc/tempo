FROM cmopipeline/facets-suite-preview:0.0.1

LABEL maintainer="Anne Marie Noronha (noronhaa@mskcc.org)" \
      contributor="Yixiao Gong (gongy@mskcc.org)" \
      contributor="Nikhil Kumar (kumarn1@mskcc.org)" \
      contributor="Philip Jonsson (jonssonp@mskcc.org)" \
      version.image="0.0.1" \
      version.facets_base="0.0.1" \
      version.htstools="0.1.1" \
      version.htslib="1.5" \
      source.htstools="https://github.com/mskcc/htstools/releases/tag/snp_pileup_0.1.1" \
      source.htslib="https://github.com/samtools/htslib/releases/tag/1.5" \
      version.annotate_with_zygosity_somatic="0.1.0" \
      version.annotate_with_zygosity_germline="0.1.0" 

ENV HTSTOOLS_VERSION 0.1.1
ENV HTSLIB_VERSION 1.5

# Install htslib
RUN cd /tmp \
    && wget https://github.com/samtools/htslib/releases/download/${HTSLIB_VERSION}/htslib-${HTSLIB_VERSION}.tar.bz2 \
    && tar xvjf htslib-${HTSLIB_VERSION}.tar.bz2 \
    && cd /tmp/htslib-${HTSLIB_VERSION} \
    && ./configure \
    && make && make install

# Install htstools
RUN cd /tmp \
    && wget https://github.com/mskcc/htstools/archive/snp_pileup_${HTSTOOLS_VERSION}.tar.gz \
    && tar xvzf snp_pileup_${HTSTOOLS_VERSION}.tar.gz \
    && cd /tmp/htstools-snp_pileup_${HTSTOOLS_VERSION} \
    && cp /tmp/htslib-${HTSLIB_VERSION}/libhts.so /usr/lib \
    && cp /tmp/htslib-${HTSLIB_VERSION}/libhts.so.2 /usr/lib \
    && cd /tmp/htstools-snp_pileup_${HTSTOOLS_VERSION} \
    && g++ -std=c++11 snp-pileup.cpp -lhts -o snp-pileup \
    && g++ -std=c++11 ppflag-fixer.cpp -lhts -o ppflag-fixer \
    # move executables into bin
    && cp /tmp/htstools-snp_pileup_${HTSTOOLS_VERSION}/snp-pileup /usr/bin \
    && cp /tmp/htstools-snp_pileup_${HTSTOOLS_VERSION}/ppflag-fixer /usr/bin

# Use old facetsSuite 1.6.3 summarize_project.py
RUN cd /usr/bin \
    && wget https://raw.githubusercontent.com/mskcc/facets-suite/1.6.3/summarize_project.py

ENV PYTHONNOUSERSITE set
ENV FACETS_OVERRIDE_EXITCODE set

COPY annotate-with-zygosity-somatic.R /usr/bin
COPY annotate-with-zygosity-germline.R /usr/bin
RUN chmod 755 /usr/bin/annotate-with-zygosity-somatic.R
RUN chmod 755 /usr/bin/annotate-with-zygosity-germline.R
