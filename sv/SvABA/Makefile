export SHELL:=/bin/bash
.ONESHELL:
export NXF_VER:=20.07.1

./nextflow:
	if module avail java/jdk1.8.0_202 1&>/dev/null; then module load java/jdk1.8.0_202; fi
	curl -fsSL get.nextflow.io | bash
install: ./nextflow



# ~~~~~ Test Workflow ~~~~~ #
run:
	if module avail java/jdk1.8.0_202 1&>/dev/null; then module load java/jdk1.8.0_202; fi
	./nextflow run -resume main.nf



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





# ~~~~~ Container ~~~~~ #
# make the Docker container
GIT_NAME:=svaba
# GIT_TAG:=$(shell git describe --tags --abbrev=0)
GIT_TAG:=dev
DOCKER_TAG:=mskcc/$(GIT_NAME):$(GIT_TAG)
docker-build:
	docker build -t "$(DOCKER_TAG)" .

# shell into the container to check that it looks right
docker-bash:
	docker run --rm -ti "$(DOCKER_TAG)" bash

# push the container to Dockerhub
# $ docker login --username=<username>
docker-push:
	docker push "$(DOCKER_TAG)"

# pull the Dockerhub container and convert to Singularity container
# NOTE: you cannot use a filename with a ':' as a Makefile target
SINGULARITY_SIF:=mskcc_$(GIT_NAME):$(GIT_TAG).sif
singularity-pull:
	unset SINGULARITY_CACHEDIR && \
	module load singularity/3.3.0 && \
	singularity pull --force --name "$(SINGULARITY_SIF)" docker://$(DOCKER_TAG)


# ~~~~~~ Reference Containers ~~~~~~~ #
# containers from the publications

# SvABA
# https://hub.docker.com/r/weischenfeldt/svaba_workflow/tags
# ea066dafc7d7
weischenfeldt.svaba_workflow.1.0.0.tar.gz:
	docker pull weischenfeldt/svaba_workflow:1.0.0 && \
	docker save weischenfeldt/svaba_workflow:1.0.0 | \
	gzip > weischenfeldt.svaba_workflow.1.0.0.tar.gz

# Delly Workflow
# https://hub.docker.com/r/weischenfeldt/bric-embl_delly_workflow/tags
# a1481f64aea7
weischenfeldt.bric-embl_delly_workflow.2.0.2_ppcg.tar.gz:
	docker pull weischenfeldt/bric-embl_delly_workflow:2.0.2_ppcg && \
	docker save weischenfeldt/bric-embl_delly_workflow:2.0.2_ppcg | \
	gzip > weischenfeldt.bric-embl_delly_workflow.2.0.2_ppcg.tar.gz

# BRASS (Sanger pipeline)
# https://hub.docker.com/r/sevenbridges/pcawg_sanger_sbg_modified
# da98ea531af7
sevenbridges.pcawg_sanger_sbg_modified.1.0.tar.gz:
	docker pull sevenbridges/pcawg_sanger_sbg_modified:1.0 && \
	docker save sevenbridges/pcawg_sanger_sbg_modified:1.0 | \
	gzip > sevenbridges.pcawg_sanger_sbg_modified.1.0.tar.gz


# PCAWG6 workflow to merge SVs from Sanger (BRASS), Broad Inst (DRanger & SnowMan) and EMBL/DKFZ (DELLY)
# https://hub.docker.com/r/weischenfeldt/pcawg_sv_merge
# ee54fa48c27b
weischenfeldt.pcawg_sv_merge.1.0.2.tar.gz:
	docker pull weischenfeldt/pcawg_sv_merge:1.0.2 && \
	docker save weischenfeldt/pcawg_sv_merge:1.0.2 | \
	gzip > weischenfeldt.pcawg_sv_merge.1.0.2.tar.gz

ALL_REF_CONTAINERS:=weischenfeldt.svaba_workflow.1.0.0.tar.gz weischenfeldt.bric-embl_delly_workflow.2.0.2_ppcg.tar.gz sevenbridges.pcawg_sanger_sbg_modified.1.0.tar.gz weischenfeldt.pcawg_sv_merge.1.0.2.tar.gz

docker-containers: $(ALL_REF_CONTAINERS)



# Singularity containers
weischenfeldt.svaba_workflow.1.0.0.sif:
	unset SINGULARITY_CACHEDIR && \
	module load singularity/3.3.0 && \
	singularity pull --force \
	--name weischenfeldt.svaba_workflow.1.0.0.sif \
	docker://weischenfeldt/svaba_workflow:1.0.0

weischenfeldt.bric-embl_delly_workflow.2.0.2_ppcg.sif:
	unset SINGULARITY_CACHEDIR && \
	module load singularity/3.3.0 && \
	singularity pull --force \
	--name weischenfeldt.bric-embl_delly_workflow.2.0.2_ppcg.sif \
	docker://weischenfeldt/bric-embl_delly_workflow:2.0.2_ppcg

sevenbridges.pcawg_sanger_sbg_modified.1.0.sif:
	unset SINGULARITY_CACHEDIR && \
	module load singularity/3.3.0 && \
	singularity pull --force \
	--name sevenbridges.pcawg_sanger_sbg_modified.1.0.sif \
	docker://sevenbridges/pcawg_sanger_sbg_modified:1.0

weischenfeldt.pcawg_sv_merge.1.0.2.sif:
	unset SINGULARITY_CACHEDIR && \
	module load singularity/3.3.0 && \
	singularity pull --force \
	--name weischenfeldt.pcawg_sv_merge.1.0.2.sif \
	docker://weischenfeldt/pcawg_sv_merge:1.0.2

ALL_SIF_CONTAINERS:=weischenfeldt.svaba_workflow.1.0.0.sif weischenfeldt.bric-embl_delly_workflow.2.0.2_ppcg.sif sevenbridges.pcawg_sanger_sbg_modified.1.0.sif weischenfeldt.pcawg_sv_merge.1.0.2.sif

singularity-containers: $(ALL_SIF_CONTAINERS)
