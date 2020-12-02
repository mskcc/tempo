export SHELL:=/bin/bash
.ONESHELL:
UNAME:=$(shell uname)
# export NXF_VER:=20.07.1
export NXF_VER:=20.01.0
export PATH:=$(CURDIR):$(CURDIR)/conda/bin:$(PATH)
unexport PYTHONPATH
unexport PYTHONHOME

./nextflow:
	if module avail java/jdk1.8.0_202 1&>/dev/null; then module load java/jdk1.8.0_202; fi
	curl -fsSL get.nextflow.io | bash


# install versions of conda for Mac or Linux, Python 2 or 3
ifeq ($(UNAME), Darwin)
CONDASH:=Miniconda3-4.7.12.1-MacOSX-x86_64.sh
endif

ifeq ($(UNAME), Linux)
CONDASH:=Miniconda3-4.7.12.1-Linux-x86_64.sh
endif

CONDAURL:=https://repo.continuum.io/miniconda/$(CONDASH)
conda:
	@echo ">>> Setting up conda..."
	@wget "$(CONDAURL)" && \
	bash "$(CONDASH)" -b -p conda && \
	rm -f "$(CONDASH)"

install: ./nextflow conda
	conda config --add channels r
	conda config --add channels bioconda
	conda install -y pysam

test-data:
	wget -O test-data.tar.gz --no-check-certificate 'https://docs.google.com/uc?export=download&id=1Jgsj4FEnwhNybaNuovyIdSXG6-iBzAIf' && \
	tar -xzvf test-data.tar.gz

# python/3.7.1
test: install test-data
	module load singularity/3.3.0 java/jdk1.8.0_202
	python3 test_pipelines.py

bash:
	module load singularity/3.3.0 java/jdk1.8.0_202
	bash

clean:
	rm -f .nextflow.log.*
	rm -f *.html.*
	rm -f trace.txt.*

clean-all: clean
	rm -f .nextflow.log*
	rm -f *.html
	rm -f trace.txt*
	rm -rf work
	rm -rf output
