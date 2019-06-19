FROM alpine:3.8

LABEL maintainer="Christopher Allan Bolipata (bolipatc@mskcc.org)" \
      contributor="Nikhil Kumar (kumarn1@mskcc.org)" \
      version.image="1.1.0" \
      version.vcf2maf="1.6.17" \
      version.vep="88" \
      version.htslib="1.9" \
      version.bcftools="1.9" \
      version.samtools="1.9" \
      version.perl="5.26.2-r1" \
      version.alpine="3.8" \
      version.oncokb_annotator="1.1.0" \
      version.filter_somatic_maf="0.2.0" \
      version.filter_germline_maf="0.1.0" \
      source.vcf2maf="https://github.com/mskcc/vcf2maf/releases/tag/v1.6.17" \
      source.vep="http://dec2016.archive.ensembl.org/info/docs/tools/vep/script/vep_download.html#versions" \
      source.htslib="https://github.com/samtools/htslib/releases/tag/1.9" \
      source.bcftools="https://github.com/samtools/bcftools/releases/tag/1.9" \
      source.samtools="https://github.com/samtools/samtools/releases/tag/1.9" \
      source.oncokb_annotator="https://github.com/oncokb/oncokb-annotator/releases/tag/v1.1.0"

ENV VCF2MAF_VERSION 1.6.17
ENV VEP_VERSION 88
ENV HTSLIB_VERSION 1.9
ENV SAMTOOLS_VERSION 1.9
ENV BCFTOOLS_VERSION 1.9
ENV ONCOKB_VERSION 1.1.0

RUN apk add --update \
      # install all the build-related tools
            && apk add ca-certificates gcc g++ make git curl curl-dev wget gzip perl perl-dev musl-dev libgcrypt-dev zlib-dev bzip2-dev xz-dev ncurses-dev python python-dev py-pip build-base freetype-dev libpng-dev openblas-dev \
      # install system packages and Perl modules
            && apk add expat-dev libressl-dev perl-net-ssleay mariadb-dev libxml2-dev perl-dbd-mysql perl-module-metadata perl-gd perl-db_file perl-archive-zip perl-cgi perl-dbi perl-encode perl-time-hires perl-file-copy-recursive perl-json \
      # install cpanminus
            && curl -L https://cpanmin.us | perl - App::cpanminus \
      # install perl libraries that VEP will need
            && cpanm --notest LWP LWP::Simple LWP::Protocol::https Archive::Extract Archive::Tar Archive::Zip \
            CGI DBI Encode version Time::HiRes File::Copy::Recursive Perl::OSType Module::Metadata \
            Sereal JSON Bio::Root::Version Set::IntervalTree PerlIO::gzip \
      # install htslib (for vep)
            && cd /tmp && wget https://github.com/samtools/htslib/releases/download/${HTSLIB_VERSION}/htslib-${HTSLIB_VERSION}.tar.bz2 \
            && tar xvjf htslib-${HTSLIB_VERSION}.tar.bz2 \
            && cd /tmp/htslib-${HTSLIB_VERSION} \
            && ./configure \
            && make && make install \
      # download/unzip vep
            && cd /tmp && wget https://github.com/Ensembl/ensembl-tools/archive/release/${VEP_VERSION}.zip \
            && unzip ${VEP_VERSION} \
      # install vep
            && cd /tmp/ensembl-tools-release-${VEP_VERSION}/scripts/variant_effect_predictor \
            && perl INSTALL.pl --AUTO a 2>&1 | tee install.log \
            && cd /tmp && mv /tmp/ensembl-tools-release-${VEP_VERSION}/scripts/variant_effect_predictor /usr/bin/vep \
      # install bcftools
            && cd /tmp && wget https://github.com/samtools/bcftools/releases/download/${BCFTOOLS_VERSION}/bcftools-${BCFTOOLS_VERSION}.tar.bz2 \
            && tar xvjf bcftools-${BCFTOOLS_VERSION}.tar.bz2 \
            && cd /tmp/bcftools-${BCFTOOLS_VERSION} \
            && make HTSDIR=/tmp/htslib-${HTSLIB_VERSION} && make install \
      # install samtools
            && cd /tmp && wget https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 \
            && tar xvjf samtools-${SAMTOOLS_VERSION}.tar.bz2 \
            && cd /tmp/samtools-${SAMTOOLS_VERSION} \
            && ./configure --with-htslib=/tmp/htslib-${HTSLIB_VERSION} \
            && make && make install \
      # install vcf2maf
            && cd /tmp && wget -O vcf2maf-v${VCF2MAF_VERSION} https://github.com/mskcc/vcf2maf/archive/v${VCF2MAF_VERSION}.zip \
            && unzip vcf2maf-v${VCF2MAF_VERSION} \
            && mkdir -p /usr/bin/vcf2maf \
            && cp -r vcf2maf-${VCF2MAF_VERSION}/* /usr/bin/vcf2maf/ \
      # install oncokb-annotator
            && cd /tmp && wget -O oncokb_annotator-v${ONCOKB_VERSION} https://github.com/oncokb/oncokb-annotator/archive/v${ONCOKB_VERSION}.zip \
            && unzip oncokb_annotator-v${ONCOKB_VERSION} \
            && mkdir -p /usr/bin/oncokb_annotator \
            && cp -r oncokb-annotator-${ONCOKB_VERSION}/* /usr/bin/oncokb_annotator/ \
      # clean up
            && rm -rf /var/cache/apk/* /tmp/*

# Install python dependencies
RUN pip install --upgrade pip
RUN pip install matplotlib

# Install R with depencies
RUN apk add --update && \
    apk add --no-cache R R-dev && \
    R -e "install.packages(c('data.table', 'rlang', 'devtools'), dependencies = TRUE, repos = 'http://cran.rstudio.com'); devtools::install_github('taylor-lab/annotateMaf')" 

# Add filter script
COPY filter-somatic-maf.R /usr/bin
COPY filter-germline-maf.R /usr/bin
RUN chmod +x /usr/bin/filter-somatic-maf.R
RUN chmod +x /usr/bin/filter-germline-maf.R

# install bash
RUN apk add bash
