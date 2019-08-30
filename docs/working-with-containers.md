# Working With Containers

Tempo relies upon containers for reproducibility and portability. In Nextflow, each process can be executed within the environment of a specified container built on the [Docker](https://www.docker.com/) or [Singularity](https://singularity.lbl.gov) platform. This page contains a basic introduction to containers and how-to's for running Tempo with containers.

* __Containers__: A container is a lightweight encapsulation of a environment to run a specific tool. One can think of this as "packaging" a tool with its dependencies into a standardized unit. The container includes code and all dependencies so the application runs quickly and reliably from one computing environment to another (e.g. system libraries, environment settings, etc.).

* __Images__: A Docker image is a building block for the container. One can think of the container as an image "put to use". Images can be hosted in a repository such as [Dockerhub](https://cloud.docker.com). When you pull the image from repository, you will have a local version of the image you can use. We you *run* the image and use it, you're using that specific container.

* __Dockerfile__: This is the recipe used to create the *image*.

* __Docker versus Singularity__: These are different container platforms. Singularity is compatible with images built on the Docker platform. For security concerns, the Juno compute cluster uses Singularity, see the [juno setup instructions](juno-setup.md). 

* __Versioning__: Docker images are versioned by their association with a *tag*. 

Tempo uses images built by Docker and hosted at [Dockerhub](https://cloud.docker.com/u/cmopipeline/). The Dockerfiles are on [GitHub](containers) and the container associated with each pipeline process is defined in the [container configuration](conf/containers.config).
