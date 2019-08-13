(window.webpackJsonp=window.webpackJsonp||[]).push([[7],{188:function(t,e,r){"use strict";r.r(e);var a=r(0),n=Object(a.a)({},function(){var t=this,e=t.$createElement,r=t._self._c||e;return r("ContentSlotsDistributor",{attrs:{"slot-key":t.$parent.slotKey}},[r("h1",{attrs:{id:"vaporware-documentation"}},[r("a",{staticClass:"header-anchor",attrs:{href:"#vaporware-documentation","aria-hidden":"true"}},[t._v("#")]),t._v(" Vaporware documentation")]),t._v(" "),r("p",[t._v("Vaporware is a computational pipeline for processing data of paired-end whole exome (WES) and whole genome sequencing (WGS) of human cancer samples with matched normals. Its components are containerized and the pipeline runs on the "),r("a",{attrs:{href:"http://hpc.mskcc.org/",target:"_blank",rel:"noopener noreferrer"}},[t._v("Juno high-performance computing cluster"),r("OutboundLink")],1),t._v(" at Memorial Sloan Kettering Cancer Cencer and on "),r("a",{attrs:{href:"https://aws.amazon.com",target:"_blank",rel:"noopener noreferrer"}},[t._v("Amazon Web Services (AWS)"),r("OutboundLink")],1),t._v(". The pipeline was written by members of the Center for Molecular Oncology.")]),t._v(" "),r("p",[t._v("These pages contain instructions on how to run the Vaporware pipeline. It also contains documentation on the bioinformatic components in the pipeline, some motivation for various parameter choices, plus sources and processing of all reference materials use.")]),t._v(" "),r("p",[t._v("If there are any questions or comments, you are welcome to "),r("a",{attrs:{href:"https://github.com/mskcc/vaporware/issues/new?title=%5BUser%20question%5D",target:"_blank",rel:"noopener noreferrer"}},[t._v("raise an issue"),r("OutboundLink")],1),t._v(".")]),t._v(" "),r("p",[r("small",[t._v("Note: Vaporware currently only supports human samples, has only been tested for exome and genome sequencing experiments, and all references files are in build GRCh37 of the human genome.")])]),t._v(" "),r("hr"),t._v(" "),r("h2",{attrs:{id:"table-of-contents"}},[r("a",{staticClass:"header-anchor",attrs:{href:"#table-of-contents","aria-hidden":"true"}},[t._v("#")]),t._v(" Table of contents")]),t._v(" "),r("h3",{attrs:{id:"_1-getting-started"}},[r("a",{staticClass:"header-anchor",attrs:{href:"#_1-getting-started","aria-hidden":"true"}},[t._v("#")]),t._v(" 1. Getting started")]),t._v(" "),r("h4",{attrs:{id:"_1-1-installation"}},[r("a",{staticClass:"header-anchor",attrs:{href:"#_1-1-installation","aria-hidden":"true"}},[t._v("#")]),t._v(" 1.1. "),r("router-link",{attrs:{to:"/installation.html"}},[t._v("Installation")])],1),t._v(" "),r("ul",[r("li",[r("router-link",{attrs:{to:"/juno-setup.html"}},[t._v("Set up on Juno")])],1),t._v(" "),r("li",[r("router-link",{attrs:{to:"/aws-setup.html"}},[t._v("Set up on AWS")])],1)]),t._v(" "),r("h4",{attrs:{id:"_1-2-usage"}},[r("a",{staticClass:"header-anchor",attrs:{href:"#_1-2-usage","aria-hidden":"true"}},[t._v("#")]),t._v(" 1.2. "),r("router-link",{attrs:{to:"/usage.html"}},[t._v("Usage")])],1),t._v(" "),r("ul",[r("li",[r("router-link",{attrs:{to:"/nextflow-basics.html"}},[t._v("Nextflow basics")])],1),t._v(" "),r("li",[r("router-link",{attrs:{to:"/working-with-containers.html"}},[t._v("Working with containers")])],1),t._v(" "),r("li",[r("router-link",{attrs:{to:"/run-pipeline.html"}},[t._v("Run pipeline")])],1)]),t._v(" "),r("h3",{attrs:{id:"_2-pipeline-contents"}},[r("a",{staticClass:"header-anchor",attrs:{href:"#_2-pipeline-contents","aria-hidden":"true"}},[t._v("#")]),t._v(" 2. Pipeline contents")]),t._v(" "),r("h4",{attrs:{id:"_2-1-bioinformatic-components"}},[r("a",{staticClass:"header-anchor",attrs:{href:"#_2-1-bioinformatic-components","aria-hidden":"true"}},[t._v("#")]),t._v(" 2.1. "),r("router-link",{attrs:{to:"/bioinformatic-components.html"}},[t._v("Bioinformatic components")])],1),t._v(" "),r("ul",[r("li",[r("router-link",{attrs:{to:"/bioinformatic-components.html#read-alignment"}},[t._v("Read alignment")])],1),t._v(" "),r("li",[r("router-link",{attrs:{to:"/bioinformatic-components.html#somatic-analyses"}},[t._v("Somatic analyses")])],1),t._v(" "),r("li",[r("router-link",{attrs:{to:"/bioinformatic-components.html#germline-analyses"}},[t._v("Germline analyses")])],1)]),t._v(" "),r("h4",{attrs:{id:"_2-2-reference-files"}},[r("a",{staticClass:"header-anchor",attrs:{href:"#_2-2-reference-files","aria-hidden":"true"}},[t._v("#")]),t._v(" 2.2. "),r("router-link",{attrs:{to:"/reference-files.html"}},[t._v("Reference files")])],1),t._v(" "),r("ul",[r("li",[r("router-link",{attrs:{to:"/intervals.html"}},[t._v("Intervals")])],1)]),t._v(" "),r("h4",{attrs:{id:"_2-3-variant-annotation-and-filtering"}},[r("a",{staticClass:"header-anchor",attrs:{href:"#_2-3-variant-annotation-and-filtering","aria-hidden":"true"}},[t._v("#")]),t._v(" 2.3. "),r("router-link",{attrs:{to:"/variant-annotation-and-filtering.html"}},[t._v("Variant annotation and filtering")])],1),t._v(" "),r("ul",[r("li",[r("router-link",{attrs:{to:"/wes-panel-of-normals.html"}},[t._v("Panel of normals for exomes")])],1)]),t._v(" "),r("h3",{attrs:{id:"_3-help-and-other-resources"}},[r("a",{staticClass:"header-anchor",attrs:{href:"#_3-help-and-other-resources","aria-hidden":"true"}},[t._v("#")]),t._v(" 3. Help and other resources")]),t._v(" "),r("ul",[r("li",[r("router-link",{attrs:{to:"/troubleshooting.html"}},[t._v("Troubleshooting")])],1),t._v(" "),r("li",[r("router-link",{attrs:{to:"/aws-glossary.html"}},[t._v("AWS glossary")])],1)]),t._v(" "),r("h3",{attrs:{id:"_4-contributing"}},[r("a",{staticClass:"header-anchor",attrs:{href:"#_4-contributing","aria-hidden":"true"}},[t._v("#")]),t._v(" 4. Contributing")]),t._v(" "),r("ul",[r("li",[r("router-link",{attrs:{to:"/contributing-to-vaporware.html"}},[t._v("Contributing to Vaporware")])],1)]),t._v(" "),r("hr")])},[],!1,null,null,null);e.default=n.exports}}]);