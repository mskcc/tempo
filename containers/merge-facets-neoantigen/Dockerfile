FROM r-base:3.5.3

RUN apt-get -y update && \
    apt-get -y install curl \ 
    openssl \
    libcurl4-openssl-dev \
    libssl-dev \
    procps


RUN R -e "install.packages('devtools')" \
    && R -e "devtools::install_github(repo = 'PoisonAlien/maftools')"

RUN R -e "install.packages('argparse')"

RUN apt-get -y install python
