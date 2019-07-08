FROM r-base:3.5.3

RUN apt-get -y update && \
    apt-get -y install curl \ 
    openssl \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    python \
    procps

RUN R -e "install.packages(c('argparse','devtools'), dependencies=TRUE)" \
    && R -e "devtools::install_github(repo = 'PoisonAlien/maftools')"

ADD merge_mafs.R /usr/local/bin
