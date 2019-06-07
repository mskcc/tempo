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
RUN tar -xvf /tmp/netMHC/netMHC-4.0a.Linux.tar.gz -C /opt \
    && sed -i 's|'/usr/cbs/packages/netMHC/4.0/netMHC-4.0'|'/opt/netMHC-4.0'|' /opt/netMHC-4.0/netMHC

# Install netMHC data
RUN cd /opt/netMHC-4.0/ \
    && wget http://www.cbs.dtu.dk/services/NetMHC-4.0/data.tar.gz \
    && gunzip -c data.tar.gz | tar xvf - \
    && rm data.tar.gz

# Install netMHC which requires modifying entries within the script
RUN tar -xvf /tmp/netMHC/netMHCpan-4.0a.Linux.tar.gz -C /opt \ 
    && sed -i 's|'/usr/cbs/packages/netMHCpan/4.0/netMHCpan-4.0'|'/opt/netMHCpan-4.0'|' /opt/netMHCpan-4.0/netMHCpan

# Install netMHCpan data
RUN cd /opt/netMHCpan-4.0/ \
    && wget http://www.cbs.dtu.dk/services/NetMHCpan-4.0/data.Linux.tar.gz \
    && gunzip -c data.Linux.tar.gz | tar xvf - \
    && rm data.Linux.tar.gz


