FROM ubuntu:18.04 

ENV TMPDIR="/tmp"

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
# NOTE: You will have to acquire them through their website, as it is not readily available for download http://www.cbs.dtu.dk/cgi-bin/nph-sw_request?netMHCpan
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
    && sed -i 's|'/usr/cbs/packages/netMHCpan/4.0/netMHCpan-4.0'|'/usr/local/bin/netMHCpan-4.0'|' /usr/local/bin/netMHCpan-4.0/netMHCpan

# Install netMHCpan data; had to make local changes (directory needed permission changes), so I'm copying them over
# This is from http://www.cbs.dtu.dk/services/NetMHCpan-4.0/data.Linux.tar.gz
COPY data.netmhcpan.tar.gz /usr/local/bin/netMHCpan-4.0/

RUN cd /usr/local/bin/netMHCpan-4.0/ \
    && gunzip -c data.netmhcpan.tar.gz | tar xvf - \
    && rm data.netmhcpan.tar.gz

# Add neoantigen-dev
COPY neoantigen-docker.config /tmp
COPY neoantigen-dev-0.0.2.tar.gz /tmp

RUN tar -xvf /tmp/neoantigen-dev-0.0.2.tar.gz -C /usr/local/bin \
    && mv /usr/local/bin/neoantigen-dev-0.0.2/ /usr/local/bin/neoantigen/ \
    && cd /usr/local/bin/neoantigen/ \
    && pip install -r requirements.txt \
    && cp /tmp/neoantigen-docker.config .
