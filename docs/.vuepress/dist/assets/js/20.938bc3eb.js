(window.webpackJsonp=window.webpackJsonp||[]).push([[20],{203:function(e,t,n){"use strict";n.r(t);var o=n(0),i=Object(o.a)({},function(){var e=this,t=e.$createElement,n=e._self._c||t;return n("ContentSlotsDistributor",{attrs:{"slot-key":e.$parent.slotKey}},[n("h1",{attrs:{id:"working-with-containers"}},[n("a",{staticClass:"header-anchor",attrs:{href:"#working-with-containers","aria-hidden":"true"}},[e._v("#")]),e._v(" Working With Containers")]),e._v(" "),n("p",[e._v("Tempo relies upon containers for reproducibility and portability. In Nextflow, each process can be executed within the environment of a specified container built on the "),n("a",{attrs:{href:"https://www.docker.com/",target:"_blank",rel:"noopener noreferrer"}},[e._v("Docker"),n("OutboundLink")],1),e._v(" or "),n("a",{attrs:{href:"https://singularity.lbl.gov",target:"_blank",rel:"noopener noreferrer"}},[e._v("Singularity"),n("OutboundLink")],1),e._v(" platform. This page contains a basic introduction to containers and how-to's for running Tempo with containers.")]),e._v(" "),n("ul",[n("li",[n("p",[n("strong",[e._v("Containers")]),e._v(': A container is a lightweight encapsulation of a environment to run a specific tool. One can think of this as "packaging" a tool with its dependencies into a standardized unit. The container includes code and all dependencies so the application runs quickly and reliably from one computing environment to another (e.g. system libraries, environment settings, etc.).')])]),e._v(" "),n("li",[n("p",[n("strong",[e._v("Images")]),e._v(': A Docker image is a building block for the container. One can think of the container as an image "put to use". Images can be hosted in a repository such as '),n("a",{attrs:{href:"https://cloud.docker.com",target:"_blank",rel:"noopener noreferrer"}},[e._v("Dockerhub"),n("OutboundLink")],1),e._v(". When you pull the image from repository, you will have a local version of the image you can use. We you "),n("em",[e._v("run")]),e._v(" the image and use it, you're using that specific container.")])]),e._v(" "),n("li",[n("p",[n("strong",[e._v("Dockerfile")]),e._v(": This is the recipe used to create the "),n("em",[e._v("image")]),e._v(".")])]),e._v(" "),n("li",[n("p",[n("strong",[e._v("Docker versus Singularity")]),e._v(": These are different container platforms. Singularity is compatible with images built on the Docker platform. For security concerns, the Juno compute cluster uses Singularity, see the "),n("router-link",{attrs:{to:"/juno-setup.html"}},[e._v("juno setup instructions")]),e._v(".")],1)]),e._v(" "),n("li",[n("p",[n("strong",[e._v("Versioning")]),e._v(": Docker images are versioned by their association with a "),n("em",[e._v("tag")]),e._v(".")])])]),e._v(" "),n("p",[e._v("Tempo uses images built by Docker and hosted at "),n("a",{attrs:{href:"https://cloud.docker.com/u/cmopipeline/",target:"_blank",rel:"noopener noreferrer"}},[e._v("Dockerhub"),n("OutboundLink")],1),e._v(". The Dockerfiles are on "),n("a",{attrs:{href:"containers"}},[e._v("GitHub")]),e._v(" and the container associated with each pipeline process is defined in the "),n("a",{attrs:{href:"conf/containers.config"}},[e._v("container configuration")]),e._v(".")])])},[],!1,null,null,null);t.default=i.exports}}]);