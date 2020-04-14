FROM rocker/tidyverse:3.6.1

LABEL maintainer="Yixiao Gong (gongy@mskcc.org)" \
      contributor="Nikhil Kumar (kumarn1@mskcc.org)" \
      contributor="Philip Jonsson (jonssonp@mskcc.org)" \
      version.image="3.0.1" \
      version.facets_suite="2.0.4" \
      version.facets="0.5.14"\
      version.alpine="3.8" \
      version.pctGCdata="0.2.0" \
      source.facets_suite="https://github.com/mskcc/facets-suite/releases/tag/2.0.4" \
      source.facets="https://github.com/mskcc/facets/archive/v0.5.14.tar.gz" \
      version.htstools="0.1.1" \
      version.htslib="1.5" \
      source.htstools="https://github.com/mskcc/htstools/releases/tag/snp_pileup_0.1.1" \
      source.htslib="https://github.com/samtools/htslib/releases/tag/1.5" \
      version.annotate_with_zygosity_somatic="0.1.0" \
      version.annotate_with_zygosity_germline="0.1.0" 

ENV HTSTOOLS_VERSION 0.1.1
ENV HTSLIB_VERSION 1.5
ENV FACETS_SUITE_VERSION 2.0.4
ENV FACETS_VERSION 0.5.14
ENV PCTGCDATA 0.2.0

# Requirements
RUN apt-get update && apt-get install -y g++ tar bzip2 libbz2-dev libc6-dev libxt-dev liblzma-dev libgtk2.0-dev libcairo2-dev xvfb xauth xfonts-base

RUN R -e "install.packages(c('Cairo','argparse','gridExtra', 'binom', 'BiocManager', 'diptest', 'egg'), repos='http://cran.us.r-project.org')" \
    && R -e "BiocManager::install('rtracklayer')"

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

# Install FACETS, pctGCdata and facets-suite
RUN cd /tmp \
	&& wget https://github.com/mskcc/facets-suite/archive/${FACETS_SUITE_VERSION}.tar.gz -O facets-suite-${FACETS_SUITE_VERSION}.tar.gz \
    && wget https://github.com/mskcc/facets/archive/v${FACETS_VERSION}.tar.gz -O facets-v${FACETS_VERSION}.tar.gz \
    && wget https://github.com/mskcc/pctGCdata/archive/v${PCTGCDATA}.tar.gz \
    && tar xvzf facets-v${FACETS_VERSION}.tar.gz \
    && tar xvzf v${PCTGCDATA}.tar.gz \
    && tar xvzf facets-suite-${FACETS_SUITE_VERSION}.tar.gz \
    && cd /tmp/pctGCdata-${PCTGCDATA} \
    && R CMD INSTALL . \
    && cd /tmp/facets-${FACETS_VERSION} \
    && R CMD INSTALL . \
    && cd /tmp/facets-suite-${FACETS_SUITE_VERSION} \
    && R CMD INSTALL . \
    && wget https://raw.githubusercontent.com/mskcc/facets-suite/1.6.3/summarize_project.py \
    && mv summarize_project.py /usr/bin \
    # correct shebang line
    && sed -i "s/opt\/common\/CentOS_6-dev\/R\/R-3.2.2\//usr\//g" *.R \
    # copy execs to /usr/bin/facets-suite
    && mkdir -p /usr/bin/facets-suite/ \
    && cp -r /tmp/facets-suite-${FACETS_SUITE_VERSION}/* /usr/bin/facets-suite/ \
    # clean up
    && rm -rf /var/cache/apk/* /tmp/*

# Use geneLevel.R old facetsSuite 1.6.3 to fix incorrect AMP filtering in geneLevel calls
RUN cd /tmp \
    && wget https://raw.githubusercontent.com/mskcc/facets-suite/1.6.3/geneLevel.R \
    && mv geneLevel.R /usr/bin/facets-suite/

ENV PYTHONNOUSERSITE set
ENV FACETS_OVERRIDE_EXITCODE set

COPY annotate-with-zygosity-somatic.R /usr/bin
COPY annotate-with-zygosity-germline.R /usr/bin
RUN chmod +x /usr/bin/annotate-with-zygosity-somatic.R
RUN chmod +x /usr/bin/annotate-with-zygosity-germline.R
