FROM ubuntu:18.04

LABEL maintainer="Allan Bolipata <bolipatc@mskcc.org>, Evan Biederstedt <evan.biederstedt@gmail.com>, Yixiao Gong <gongy@mskcc.org>" \
    version.image="0.3.2-hotfix" \
    version.neoantigen-dev="0.3.2-hotfix" \
    version.netMHC="4.0a" \
    version.netMHCpan="4.0a" \
    version.python="2.7.15rc1"

ENV TMPDIR="/tmp"
ENV NEOANTIGEN_VERSION 0.3.2-hotfix

RUN apt-get update && apt-get install -y \
    tcsh \
    libnss-sss \
    git \
    python-pip \
    build-essential \
    wget \
    &&  apt-get clean && apt-get purge \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# For copying netMHC and netMHCpan gz files for installation. This Dockerfile expects them to be in the same directory
# NOTE: You will have to acquire them through their website, as they are not readily available for download from:
#https://services.healthtech.dtu.dk/software.php
COPY netMHC*.gz /tmp/netMHC/

# Install netMHC, which requires modifying entries within the script
RUN tar -xvf /tmp/netMHC/netMHC-4.0a.Linux.tar.gz -C /usr/local/bin \
    && sed -i 's|'/usr/cbs/packages/netMHC/4.0/netMHC-4.0'|'/usr/local/bin/netMHC-4.0'|' /usr/local/bin/netMHC-4.0/netMHC

# Install netMHC data
RUN cd /usr/local/bin/netMHC-4.0/ \
    && wget http://www.cbs.dtu.dk/services/NetMHC-4.0/data.tar.gz \
    && gunzip -c data.tar.gz | tar xvf - \
    && rm data.tar.gz

# Install netMHCpan, which requires modifying entries within the script; also make sure the permissions for the stuff inside the tar.gz is at least drwxr-xr-x
RUN tar -xvf /tmp/netMHC/netMHCpan-4.0a.Linux.tar.gz -C /usr/local/bin \
    && sed -i 's|'/usr/cbs/packages/netMHCpan/4.0/netMHCpan-4.0'|'/usr/local/bin/netMHCpan-4.0'|' /usr/local/bin/netMHCpan-4.0/netMHCpan \
    && chmod -R o+rx /usr/local/bin/netMHCpan-4.0

# Install netMHCpan data; directory needed permission changes, drwxr-xr-x
RUN cd /usr/local/bin/netMHCpan-4.0/ \
    && wget http://www.cbs.dtu.dk/services/NetMHCpan-4.0/data.Linux.tar.gz \
    && tar -zxvf data.Linux.tar.gz \
    && chmod -R o+rx data/ \
    && rm -rf data.Linux.tar.gz

# Add neoantigen-dev
RUN cd /tmp \
    && wget https://github.com/taylor-lab/neoantigen-dev/archive/${NEOANTIGEN_VERSION}.tar.gz -O neoantigen-dev-${NEOANTIGEN_VERSION}.tar.gz

RUN tar -xvf /tmp/neoantigen-dev-${NEOANTIGEN_VERSION}.tar.gz -C /usr/local/bin \
    && mv /usr/local/bin/neoantigen-dev-${NEOANTIGEN_VERSION}/ /usr/local/bin/neoantigen/ \
    && cd /usr/local/bin/neoantigen/ \
    && pip install -r requirements.txt
