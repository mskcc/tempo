FROM rocker/tidyverse:4.1.3


LABEL maintainer="Anne Marie Noronha (noronhaa@mskcc.org)" \
    version.base="4.1.3" \
    version.image="0.0.2"

RUN apt-get update && apt-get install -y \
	zlib1g-dev \
	libbz2-dev \
	liblzma-dev \
	libc6-dev \
	build-essential \
	git

RUN R -e "install.packages(c('roxygen2','getopt'))"
RUN R -e "install.packages('devtools')"
RUN R -e "devtools::install_github('r-lib/xml2@v1.3.3')"
RUN R -e "install.packages(c('XML','rversions','data.table'))"
RUN R -e "devtools::install_github('linxihui/NNLM@0.4.4') "
RUN R -e "install.packages('BiocManager') ; BiocManager::install(c('BSgenome.Mmusculus.UCSC.mm10','BSgenome.Hsapiens.1000genomes.hs37d5','BSgenome.Hsapiens.UCSC.hg38'))"
RUN apt-get install -y libgmp3-dev
RUN R -e "options(download.file.method = 'wget') ;devtools::install_github('Nik-Zainal-Group/signature.tools.lib@v2.1.2')"

