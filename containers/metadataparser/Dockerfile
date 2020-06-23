FROM amancevice/pandas:0.24.2

LABEL maintainer="Evan Biederstedt <evan.biederstedt@gmail.com>" \
    contributor="Philip Jonsson <philip.jonsson@mskcc.org>; Yixiao Gong <gongy@mskcc.org>; Anne Marie Noronha <noronhaa@mskcc.org>" \
    version.image="0.5.4" \
    version.create_metadata_file="0.2.3" \
    version.create-aggregate-qc-file="0.2.0" \
    description="Debian image to run script with Python3.x"

RUN apt-get update && apt-get install bedtools

RUN pip install --upgrade pip
RUN pip install 'pybedtools>=0.8.0'
RUN pip install 'joblib>=0.12.5'

RUN apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get install --yes \
    r-base-core r-base-dev
RUN R -e "install.packages(c('data.table'), repos = 'http://cran.r-project.org')"

COPY create_metadata_file.py /usr/bin
COPY create-aggregate-qc-file.R /usr/bin
RUN chmod +rx /usr/bin/create_metadata_file.py
RUN chmod +rx /usr/bin/create-aggregate-qc-file.R
