FROM rocker/r-ubuntu:18.04

LABEL maintainer="Anne Marie Noronha <noronhaa@mskcc.org>" \
      version.ubuntu="18.04" \
      version.image="0.0.1"

RUN apt-get update && apt-get install -y \
    apt-utils \
    build-essential \
    git \
    libgsl23 \
    libgsl-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libgmp3-dev

RUN R -e "install.packages('roxygen2')"
RUN R -e "install.packages('devtools')"
RUN R -e "devtools::install_github('r-lib/xml2')"
RUN R -e "install.packages(c('XML','rversions','data.table'))"
RUN R -e "devtools::install_github('linxihui/NNLM') " 
RUN R -e "install.packages('BiocManager') ; BiocManager::install(c('BSgenome.Mmusculus.UCSC.mm10','BSgenome.Hsapiens.1000genomes.hs37d5','BSgenome.Hsapiens.UCSC.hg38'))"
RUN R -e "options(download.file.method = 'wget') ;devtools::install_github('Nik-Zainal-Group/signature.tools.lib')"

