FROM amancevice/pandas:0.24.2

LABEL maintainer="Evan Biederstedt <evan.biederstedt@gmail.com>" \
    version.image="0.2.0" \
    description="Debian image to run script with Python3.x"

RUN apt-get update && apt-get install bedtools

RUN pip install --upgrade pip && \ 
    pip install 'pybedtools>=0.8.0'

ADD create_metadata_file.py /usr/bin/create_metadata_file.py


