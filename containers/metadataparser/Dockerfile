FROM amancevice/pandas:0.24.2

LABEL maintainer="Evan Biederstedt <evan.biederstedt@gmail.com>" \
    contributor="Philip Jonsson <philip.jonsson@mskcc.org>" \
    version.image="0.3.0" \
    description="Debian image to run script with Python3.x"

RUN apt-get update && apt-get install bedtools

RUN pip install --upgrade pip
RUN pip install 'pybedtools>=0.8.0'
RUN pip install 'joblib>=0.12.5'

COPY create_metadata_file.py /usr/bin
COPY create-aggregate-qc-file.py /usr/bin
RUN chmod +x /usr/bin/create_metadata_file.py
RUN chmod +x /usr/bin/create-aggregate-qc-file.py
