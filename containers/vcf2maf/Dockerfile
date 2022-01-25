FROM ubuntu:18.04

LABEL maintainer="Christopher Allan Bolipata (bolipatc@mskcc.org)" \
      contributor="Nikhil Kumar (kumarn1@mskcc.org)" \
      contributor="Philip Jonsson (jonssonp@mskcc.org)" \
      contributor="Anne Marie Noronha (noronhaa@mskcc.org)" \
      version.image="1.2.7" \
      version.vcf2maf="1.6.17" \
      version.vep="88" \
      version.htslib="1.9" \
      version.bcftools="1.9" \
      version.samtools="1.9" \
      version.perl="5.26.2-r1" \
      version.alpine="3.8" \
      version.oncokb_annotator="1.1.0" \
      version.filter_somatic_maf="0.6.2" \
      version.filter_germline_maf="0.2.2" \
      version.annotateMaf="1.0.2" \
      source.vcf2maf="https://github.com/mskcc/vcf2maf/releases/tag/v1.6.17" \
      source.vep="http://dec2016.archive.ensembl.org/info/docs/tools/vep/script/vep_download.html#versions" \
      source.htslib="https://github.com/samtools/htslib/releases/tag/1.9" \
      source.bcftools="https://github.com/samtools/bcftools/releases/tag/1.9" \
      source.samtools="https://github.com/samtools/samtools/releases/tag/1.9" \
      source.oncokb_annotator="https://github.com/oncokb/oncokb-annotator/releases/tag/v1.1.0" \
      source.annotateMaf="https://github.com/taylor-lab/annotateMaf"

ENV VCF2MAF_VERSION 1.6.17
ENV VEP_VERSION 88
ENV HTSLIB_VERSION 1.9
ENV SAMTOOLS_VERSION 1.9
ENV BCFTOOLS_VERSION 1.9
ENV ONCOKB_VERSION 1.1.0
ENV ANNOTATE_MAF 1.0.1

RUN apt-get update --fix-missing && apt-get install -y \
    tcsh \
    libnss-sss \
    git \
    python2.7 \
    python-pip \
    build-essential \
    wget \
    && apt-get clean && apt-get purge \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN apt-get update && apt-get install -y \
    cpanminus \
    perl \
    curl \
    libssl-dev \
    libperlio-gzip-perl \
    libgd-perl \
    && apt-get clean && apt-get purge \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Install python dependencies
RUN pip install --upgrade pip
RUN pip install matplotlib

RUN apt-get update && apt-get install -y libnet-ssleay-perl libcrypt-ssleay-perl libextutils-pkgconfig-perl

# install perl libraries that VEP will need
RUN cpanm --notest LWP LWP::Simple LWP::Protocol::https 

RUN cpanm --notest Archive::Extract Archive::Tar Archive::Zip

RUN cpanm --notest CGI DBI Encode version Time::HiRes File::Copy::Recursive Perl::OSType Module::Metadata

RUN cpanm --notest Sereal JSON Bio::Root::Version Set::IntervalTree 

RUN apt-get update && apt-get install -y zlib1g-dev bzip2 libbz2-dev liblzma-dev unzip

# install htslib (for vep)
RUN cd /tmp && wget https://github.com/samtools/htslib/releases/download/${HTSLIB_VERSION}/htslib-${HTSLIB_VERSION}.tar.bz2 \
            && tar xvjf htslib-${HTSLIB_VERSION}.tar.bz2 \
            && cd /tmp/htslib-${HTSLIB_VERSION} \
            && ./configure \
            && make && make install \
            && rm -rf /tmp/*

# install vep
RUN cd /tmp && wget https://github.com/Ensembl/ensembl-tools/archive/release/${VEP_VERSION}.zip \
            && unzip ${VEP_VERSION} \
            && cd /tmp/ensembl-tools-release-${VEP_VERSION}/scripts/variant_effect_predictor \
            && perl INSTALL.pl --AUTO a 2>&1 | tee install.log \
            && cd /tmp && mv /tmp/ensembl-tools-release-${VEP_VERSION}/scripts/variant_effect_predictor /usr/bin/vep \
            && rm -rf /tmp/*

# Install samtools
RUN apt-get update && \
    apt-get install --yes \
    libncurses5-dev \
    vcftools \
    libbz2-dev \
    libxml2 \
    libxml2-dev \
    libcurl4-openssl-dev \
    liblzma-dev && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

RUN cd /tmp && curl -L -o tmp.tar.gz https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 && \
    mkdir samtools && \
    tar -C samtools --strip-components 1 -jxf tmp.tar.gz && \
    cd samtools && \
    ./configure && \
    make && \
    make install && \
    cd .. && \
    rm -rf /tmp/*

RUN cd /tmp && curl -L -o tmp2.tar.gz https://github.com/samtools/bcftools/releases/download/${BCFTOOLS_VERSION}/bcftools-${BCFTOOLS_VERSION}.tar.bz2 && \
    mkdir bcftools && \
    tar -C bcftools --strip-components 1 -jxf tmp2.tar.gz && \
    cd bcftools && \
    make && \
    make install && \
    cd .. && \
    rm -rf /tmp/*

# Install vcf2maf
WORKDIR /opt/

RUN cd /opt/ && curl -ksSL -o tmp.tar.gz https://github.com/mskcc/vcf2maf/archive/v${VCF2MAF_VERSION}.tar.gz && \
    tar --strip-components 1 -zxf tmp.tar.gz && \
    rm tmp.tar.gz && \
    chmod +x *.pl

# Install R with depencies
RUN apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get install --yes \
    r-base-core r-base-dev

# versioning R packages for reproducibility and compatibility
RUN R -e "install.packages('remotes')"
COPY r_versions.txt /
RUN R -e "require(remotes) ; x <- read.table('/r_versions.txt',sep='\t',header=T) ; for (i in 1:dim(x)[[1]]){ install_version(as.character(x\$package[i]), version = as.character(x\$version[i]), repos = 'http://cran.us.r-project.org') }" 

RUN R -e "devtools::install_github('taylor-lab/annotateMaf', ref = '${ANNOTATE_MAF}', upgrade = 'never')"

# install oncokb-annotator
RUN cd /tmp && wget -O oncokb_annotator-v${ONCOKB_VERSION} https://github.com/oncokb/oncokb-annotator/archive/v${ONCOKB_VERSION}.zip \
            && unzip oncokb_annotator-v${ONCOKB_VERSION} \
            && mkdir -p /usr/bin/oncokb_annotator \
            && cp -r oncokb-annotator-${ONCOKB_VERSION}/* /usr/bin/oncokb_annotator/ \
            && rm -rf /var/cache/apk/* /tmp/*
   
# Add filter script
COPY filter-somatic-maf.R /usr/bin
COPY filter-germline-maf.R /usr/bin
RUN chmod 755 /usr/bin/filter-somatic-maf.R
RUN chmod 755 /usr/bin/filter-germline-maf.R

