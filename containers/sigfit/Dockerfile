FROM r-base:3.5.3

LABEL authors="Jason Hwee (hweej@mskcc.org)" \
      version.image="1.0.0" \
      version.sigfit="1.3.1" \
      version.r-base="3.5.3" \
      source.sigfit="https://github.com/kgori/sigfit/releases/tag/1.3.1"

RUN apt-get -y update && \
    apt-get -y install curl \ 
    openssl \
    libcurl4-openssl-dev \
    libssl-dev

# Add sigfit dependencies && install sigfit
RUN R -e "install.packages('rstan', repos = 'https://cloud.r-project.org/', dependencies = TRUE)" \
    && R -e "install.packages('devtools')" \ 
    && R -e "devtools::install_github('kgori/sigfit', args = '--preclean', build_vignettes = TRUE)"

# Generic rscripts location, can be expanded
COPY rscripts /rscripts
ENV PATH="/rscripts:${PATH}"
WORKDIR /data

# ENTRYPOINT ["R"]
# CMD ["R"]
