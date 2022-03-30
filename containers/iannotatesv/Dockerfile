FROM python:2.7.15

LABEL maintainer="Anne Marie Noronha (noronhaa@mskcc.org)" \
	version.image="0.0.2"

RUN python -m pip install \
	pandas==0.24.2 \
	biopython==1.76 \
	Pillow==6.2.1 \
	openpyxl==2.6.4 \
	reportlab==3.5.2 \
	coloredlogs==14.0

RUN apt-get install -y git 

RUN curl -s https://packagecloud.io/install/repositories/github/git-lfs/script.deb.sh | bash && apt-get install git-lfs

RUN mkdir -p /usr/bin/ && cd /usr/bin/ && git clone https://github.com/rhshah/iAnnotateSV.git \
	&& cd iAnnotateSV \
	&& git reset --hard a2f86543925169219c91fe4e3de5412a69f735a4

RUN mkdir -p /tmp \
    && apt-get update \
    && apt-get install -y \
    bedtools
RUN python -m pip install pybedtools==0.8.2
