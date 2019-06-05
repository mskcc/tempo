FROM centos:7.6.1810

LABEL maintainer="Jason Hwee <hweej@mskcc.org>" \
    version.image="0.1.0" \
    version.centos="7-6.1810" \
    version.R="3.5.2" \
    version.python="2.7.5" \
    version.msig="0.19.7" \
    ## source.msig="https://github.com/hweej/mutation-signatures.git"
    source.msig="https://github.com/mskcc/mutation-signatures.git"

RUN yum -y update \ 
    && yum -y install python2-devel \ 
    gcc \
    tkinter \
    git \
    epel-release 
RUN yum -y install BEDTools \
    R 

RUN R -e "install.packages('data.table', dependencies=TRUE, repos = 'https://Rdatatable.github.io/data.table')"

# Install Pip
RUN curl -O "https://bootstrap.pypa.io/get-pip.py" \ 
    && python get-pip.py

# Fetch Mutation-signatures package
## transferred Jason's branch to github.com/mskcc
## RUN git clone --single-branch --branch feature/make-nice-for-docker https://github.com/hweej/mutation-signatures.git

RUN git clone --single-branch --branch feature/make-nice-for-docker https://github.com/mskcc/mutation-signatures.git  

RUN pip install -r mutation-signatures/requirements.txt

ENV PATH="/mutation-signatures:${PATH}"
ENV PYTHONPATH="/mutation-signatures"

WORKDIR /data


