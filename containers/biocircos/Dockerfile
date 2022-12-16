FROM rocker/tidyverse:4.2.2


LABEL maintainer="Anne Marie Noronha (noronhaa@mskcc.org)" \
    version.base="4.2.2" \
    version.image="0.0.1"

RUN apt-get update && apt-get install -y \
	zlib1g-dev \
	libbz2-dev \
	liblzma-dev \
	libc6-dev \
	build-essential \
	git \
	ghostscript \
	libgmp3-dev

RUN R -e "install.packages(c('BioCircos','data.table'))"
RUN R -e "install.packages('getopt')"
RUN R -e "install.packages('htmlwidgets')"

