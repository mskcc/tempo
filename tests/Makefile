export SHELL:=/bin/bash
.ONESHELL:
# export NXF_VER:=20.07.1
export NXF_VER:=20.01.0
export PATH:=$(CURDIR):$(PATH)

./nextflow:
	if module avail java/jdk1.8.0_202 1&>/dev/null; then module load java/jdk1.8.0_202; fi
	curl -fsSL get.nextflow.io | bash
install: ./nextflow

test-data.tar.gz:
	wget -O test-data.tar.gz --no-check-certificate 'https://docs.google.com/uc?export=download&id=1Jgsj4FEnwhNybaNuovyIdSXG6-iBzAIf'

test-data: test-data.tar.gz
	tar -xzvf test-data.tar.gz

test:
	module load python/3.7.1 singularity/3.3.0 java/jdk1.8.0_202
	python3 test_pipelines.py

bash:
	module load python/3.7.1 singularity/3.3.0 java/jdk1.8.0_202
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
