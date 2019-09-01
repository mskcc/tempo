FROM frolvlad/alpine-glibc:alpine-3.9

LABEL maintainer="Yixiao Gong (gongy@mskcc.org)" \
      contributor="Nikhil Kumar (kumarn1@mskcc.org)" \
      contributor="Philip Jonsson (jonssonp@mskcc.org)" \
      version.image="2.2.2" \
      version.facets_suite="1.6.2" \
      version.facets="0.5.14"\
      version.alpine="3.8" \
      version.r="3.5.1" \
      version.pctGCdata="0.2.0" \
      source.facets_suite="https://github.com/mskcc/facets-suite/releases/tag/1.6.0" \
      source.facets="https://github.com/mskcc/facets/archive/v0.5.14.tar.gz" \
      source.r="https://pkgs.alpinelinux.org/package/edge/community/x86/R" \
      version.htstools="0.1.1" \
      version.htslib="1.5" \
      source.htstools="https://github.com/mskcc/htstools/releases/tag/snp_pileup_0.1.1" \
      source.htslib="https://github.com/samtools/htslib/releases/tag/1.5" \
      version.annotate_with_zygosity_somatic="0.1.0" \
      version.annotate_with_zygosity_germline="0.1.0" 

ENV HTSTOOLS_VERSION 0.1.1
ENV HTSLIB_VERSION 1.5
ENV FACETS_SUITE_VERSION 1.6.2
ENV FACETS_VERSION 0.5.14
ENV PCTGCDATA 0.2.0

RUN apk add --update \
    # for wget
        && apk add ca-certificates openssl \
    # for compilation
        && apk add bash build-base musl-dev python py-pip python-dev \
        && apk add gcc g++ make git curl curl-dev wget gzip libgcrypt-dev zlib-dev bzip2-dev xz-dev ncurses-dev argp-standalone \
    # download, unzip, install htslib
        && cd /tmp \
        && wget https://github.com/samtools/htslib/releases/download/${HTSLIB_VERSION}/htslib-${HTSLIB_VERSION}.tar.bz2 \
        && tar xvjf htslib-${HTSLIB_VERSION}.tar.bz2 \
        && cd /tmp/htslib-${HTSLIB_VERSION} \
        && ./configure \
        && make && make install \
    # download, unzip
        && cd /tmp \
        && wget https://github.com/mskcc/htstools/archive/snp_pileup_${HTSTOOLS_VERSION}.tar.gz \
        && tar xvzf snp_pileup_${HTSTOOLS_VERSION}.tar.gz \
        && cd /tmp/htstools-snp_pileup_${HTSTOOLS_VERSION} \
        && cp /tmp/htslib-${HTSLIB_VERSION}/libhts.so /usr/lib \
        && cp /tmp/htslib-${HTSLIB_VERSION}/libhts.so.2 /usr/lib \
    # install snp-pileup and ppflag-fixer
        && cd /tmp/htstools-snp_pileup_${HTSTOOLS_VERSION} \
        && g++ -std=c++11 snp-pileup.cpp -largp -lhts -o snp-pileup \
        && g++ -std=c++11 ppflag-fixer.cpp -largp -lhts -o ppflag-fixer \
    # move executables into bin
        && cp /tmp/htstools-snp_pileup_${HTSTOOLS_VERSION}/snp-pileup /usr/bin \
        && cp /tmp/htstools-snp_pileup_${HTSTOOLS_VERSION}/ppflag-fixer /usr/bin \
    # install cairo and R package system dependencies
        && apk add cairo cairo-dev libxt-dev libxml2-dev font-xfree86-type1 \
    # download and install R
        && apk add R R-dev \
    # change permissions for /tmp; still doesn't resolve R warnings
        && chmod  777 -R /tmp \
    # download and unzip facets, facets-suite, pctGCdata
        && cd /tmp \
	&& wget https://github.com/mskcc/facets-suite/archive/${FACETS_SUITE_VERSION}.tar.gz -O facets-suite-${FACETS_SUITE_VERSION}.tar.gz \
    # download Facets
        && wget https://github.com/mskcc/facets/archive/v${FACETS_VERSION}.tar.gz -O facets-v${FACETS_VERSION}.tar.gz \
    # download pctGCdata
        && wget https://github.com/mskcc/pctGCdata/archive/v${PCTGCDATA}.tar.gz \
    # unpack Facets
        && tar xvzf facets-v${FACETS_VERSION}.tar.gz \
    # unpack pctGCdata
        && tar xvzf v${PCTGCDATA}.tar.gz \
    # unpack Facets_Suite
        && tar xvzf facets-suite-${FACETS_SUITE_VERSION}.tar.gz \
    # install
    # install pctGCdata
        && cd /tmp/pctGCdata-${PCTGCDATA} \
        && R CMD INSTALL . \
    # install Facets
        && cd /tmp/facets-${FACETS_VERSION} \
        && R CMD INSTALL . \
    # install Facets-Suite R dependencies
        && R -e "install.packages(c('ggplot2','Cairo','argparse','data.table','gtable','gridExtra','bit64','plyr', 'dplyr','tidyr','jsonlite','stringr','reshape2', 'binom'), repos='http://cran.us.r-project.org')" \
        && R -e "source('https://bioconductor.org/biocLite.R'); biocLite('rtracklayer')" \ 
    # install Facets-Suite
        && cd /tmp/facets-suite-${FACETS_SUITE_VERSION} \
    # correct shebang line
        && sed -i "s/opt\/common\/CentOS_6-dev\/R\/R-3.2.2\//usr\//g" *.R \
    # copy execs to /usr/bin/facets-suite
        && mkdir -p /usr/bin/facets-suite/ \
        && cp -r /tmp/facets-suite-${FACETS_SUITE_VERSION}/* /usr/bin/facets-suite/ \
    # clean up
        && rm -rf /var/cache/apk/* /tmp/*

ENV PYTHONNOUSERSITE set
ENV FACETS_OVERRIDE_EXITCODE set

COPY annotate-with-zygosity-somatic.R /usr/bin
COPY annotate-with-zygosity-germline.R /usr/bin
RUN chmod +x /usr/bin/annotate-with-zygosity-somatic.R
RUN chmod +x /usr/bin/annotate-with-zygosity-germline.R
