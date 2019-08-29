(window.webpackJsonp=window.webpackJsonp||[]).push([[4],{184:function(e,t,r){e.exports=r.p+"assets/img/pipeline-flowchart.ca996bb1.png"},191:function(e,t,r){"use strict";r.r(t);var n=r(0),a=Object(n.a)({},function(){var e=this,t=e.$createElement,n=e._self._c||t;return n("ContentSlotsDistributor",{attrs:{"slot-key":e.$parent.slotKey}},[n("h1",{attrs:{id:"bioinformatic-components"}},[n("a",{staticClass:"header-anchor",attrs:{href:"#bioinformatic-components","aria-hidden":"true"}},[e._v("#")]),e._v(" Bioinformatic Components")]),e._v(" "),n("p",[e._v("The three main functions of the pipeline are:")]),e._v(" "),n("ol",[n("li",[e._v("Read alignment")]),e._v(" "),n("li",[e._v("Somatic variant detection")]),e._v(" "),n("li",[e._v("Germline variant detection")])]),e._v(" "),n("p",[e._v("Additionally, various QC metrics are generated. Below are described the separate modules tools used. The following diagram outlines the workflow:")]),e._v(" "),n("img",{attrs:{id:"diagram",src:r(184)}}),e._v(" "),n("p",[n("small",[e._v("Note: The pipeline can be run with already-aligned BAM files as input, which avoids the first of these three modules.")])]),e._v(" "),n("h2",{attrs:{id:"read-alignment"}},[n("a",{staticClass:"header-anchor",attrs:{href:"#read-alignment","aria-hidden":"true"}},[e._v("#")]),e._v(" Read Alignment")]),e._v(" "),n("p",[e._v("Tempo accepts as input sequencing reads from one or multiple FASTQ file pairs (corresponding to separate sequencing lanes) per sample, as "),n("router-link",{attrs:{to:"/running-the-pipeline.html#the-mapping-file"}},[e._v("described")]),e._v(". These are aligned against the human genome using common practices, which include:")],1),e._v(" "),n("ul",[n("li",[n("strong",[e._v("Alignment")]),e._v(" using "),n("a",{attrs:{href:"http://bio-bwa.sourceforge.net/",target:"_blank",rel:"noopener noreferrer"}},[e._v("BWA mem"),n("OutboundLink")],1),e._v(", followed by conversion to BAM file format and sorting using "),n("a",{attrs:{href:"https://samtools.github.io",target:"_blank",rel:"noopener noreferrer"}},[e._v("samtools"),n("OutboundLink")],1),e._v(".")]),e._v(" "),n("li",[n("strong",[e._v("Merging")]),e._v(" of BAM files across sequencing lanes using "),n("a",{attrs:{href:"https://samtools.github.io",target:"_blank",rel:"noopener noreferrer"}},[e._v("samtools"),n("OutboundLink")],1),e._v(".")]),e._v(" "),n("li",[n("strong",[e._v("PCR-duplicate marking")]),e._v(" using "),n("a",{attrs:{href:"https://software.broadinstitute.org/gatk",target:"_blank",rel:"noopener noreferrer"}},[e._v("GATK MarkDuplicates"),n("OutboundLink")],1),e._v(".")]),e._v(" "),n("li",[n("strong",[e._v("Base-quality score recalibration")]),e._v(" with "),n("a",{attrs:{href:"https://software.broadinstitute.org/gatk/",target:"_blank",rel:"noopener noreferrer"}},[e._v("GATK BaseRecalibrator and ApplyBQSR"),n("OutboundLink")],1),e._v(".")])]),e._v(" "),n("h2",{attrs:{id:"somatic-analyses"}},[n("a",{staticClass:"header-anchor",attrs:{href:"#somatic-analyses","aria-hidden":"true"}},[e._v("#")]),e._v(" Somatic Analyses")]),e._v(" "),n("ul",[n("li",[n("strong",[e._v("SNVs and indels")]),e._v(" are called using "),n("a",{attrs:{href:"https://software.broadinstitute.org/gatk/documentation/tooldocs/4.beta.4/org_broadinstitute_hellbender_tools_walkers_mutect_Mutect2.php",target:"_blank",rel:"noopener noreferrer"}},[e._v("MuTect2"),n("OutboundLink")],1),e._v(" and "),n("a",{attrs:{href:"https://github.com/Illumina/strelka",target:"_blank",rel:"noopener noreferrer"}},[e._v("Strelka2"),n("OutboundLink")],1),e._v(". Subsequently, they are combined, annotated and filtered as described "),n("router-link",{attrs:{to:"/variant-annotation-and-filtering.html#somatic-snvs-and-indels"}},[e._v("in the section on variant annotation and filtering")]),e._v(".")],1),e._v(" "),n("li",[n("strong",[e._v("Structural variants")]),e._v(" are detected by "),n("a",{attrs:{href:"https://github.com/dellytools/delly",target:"_blank",rel:"noopener noreferrer"}},[e._v("Delly"),n("OutboundLink")],1),e._v(" and "),n("a",{attrs:{href:"https://github.com/Illumina/manta",target:"_blank",rel:"noopener noreferrer"}},[e._v("Manta"),n("OutboundLink")],1),e._v(" then combined, filtered and annotated as described "),n("router-link",{attrs:{to:"/variant-annotation-and-filtering.html#somatic-and-germline-svs"}},[e._v("in the section on variant annotation and filtering")]),e._v(".")],1),e._v(" "),n("li",[n("strong",[e._v("Copy-number analysis")]),e._v(" is performed with "),n("a",{attrs:{href:"https://github.com/mskcc/facets",target:"_blank",rel:"noopener noreferrer"}},[e._v("FACETS"),n("OutboundLink")],1),e._v(" and processed using "),n("a",{attrs:{href:"https://github.com/mskcc/facets-suite",target:"_blank",rel:"noopener noreferrer"}},[e._v("facets-suite"),n("OutboundLink")],1),e._v(". Locus-specific copy-number, purity and ploidy estimates are integrated with the SNV/indel calls to perform clonality and zygosity analyses.")]),e._v(" "),n("li",[n("strong",[e._v("Microsatellite instability")]),e._v(" is detected using "),n("a",{attrs:{href:"https://github.com/ding-lab/msisensor",target:"_blank",rel:"noopener noreferrer"}},[e._v("MSIsensor"),n("OutboundLink")],1),e._v(".")]),e._v(" "),n("li",[n("strong",[e._v("HLA genotyping")]),e._v(" is performed with "),n("a",{attrs:{href:"https://software.broadinstitute.org/cancer/cga/polysolver",target:"_blank",rel:"noopener noreferrer"}},[e._v("POLYSOLVER"),n("OutboundLink")],1),e._v(".")]),e._v(" "),n("li",[n("strong",[e._v("LOH at HLA loci")]),e._v(" is assessed with "),n("a",{attrs:{href:"https://github.com/mskcc/lohhla",target:"_blank",rel:"noopener noreferrer"}},[e._v("LOHHLA"),n("OutboundLink")],1),e._v(".")]),e._v(" "),n("li",[n("strong",[e._v("Mutational signatures")]),e._v(" are inferred with "),n("a",{attrs:{href:"https://github.com/mskcc/mutation-signatures",target:"_blank",rel:"noopener noreferrer"}},[e._v("https://github.com/mskcc/mutation-signatures"),n("OutboundLink")],1),e._v(".")]),e._v(" "),n("li",[n("strong",[e._v("Neoantigen prediction")]),e._v(" using estimates of class I MHC binding affinity is performed with "),n("a",{attrs:{href:"https://www.ncbi.nlm.nih.gov/pubmed/28978689",target:"_blank",rel:"noopener noreferrer"}},[e._v("NetMHC 4.0"),n("OutboundLink")],1),e._v(" and integrated into the set of SNV/indel calls using "),n("a",{attrs:{href:"https://github.com/taylor-lab/neoantigen-dev",target:"_blank",rel:"noopener noreferrer"}},[e._v("https://github.com/taylor-lab/neoantigen-dev"),n("OutboundLink")],1),e._v(".")])]),e._v(" "),n("h2",{attrs:{id:"germline-analyses"}},[n("a",{staticClass:"header-anchor",attrs:{href:"#germline-analyses","aria-hidden":"true"}},[e._v("#")]),e._v(" Germline Analyses")]),e._v(" "),n("ul",[n("li",[n("strong",[e._v("SNVs and indels")]),e._v(" are called using "),n("a",{attrs:{href:"https://software.broadinstitute.org/gatk/documentation/tooldocs/4.0.8.0/org_broadinstitute_hellbender_tools_walkers_haplotypecaller_HaplotypeCaller.php",target:"_blank",rel:"noopener noreferrer"}},[e._v("HaplotypeCaller"),n("OutboundLink")],1),e._v(" and "),n("a",{attrs:{href:"https://github.com/Illumina/strelka",target:"_blank",rel:"noopener noreferrer"}},[e._v("Strelka2"),n("OutboundLink")],1),e._v(". Subsequently, they are combined, annotated and filtered as described "),n("router-link",{attrs:{to:"/variant-annotation-and-filtering.html#germline-snvs-and-indels"}},[e._v("in the section on variant annotation and filtering")]),e._v(".")],1),e._v(" "),n("li",[n("strong",[e._v("Structural variants")]),e._v(" are detected by "),n("a",{attrs:{href:"https://github.com/dellytools/delly",target:"_blank",rel:"noopener noreferrer"}},[e._v("Delly"),n("OutboundLink")],1),e._v(" and "),n("a",{attrs:{href:"https://github.com/Illumina/manta",target:"_blank",rel:"noopener noreferrer"}},[e._v("Manta"),n("OutboundLink")],1),e._v(" then combined, filtered and annotated as described "),n("router-link",{attrs:{to:"/variant-annotation-and-filtering.html#somatic-and-germline-svs"}},[e._v("in the section on variant annotation and filtering")]),e._v(".")],1)]),e._v(" "),n("h2",{attrs:{id:"quality-control"}},[n("a",{staticClass:"header-anchor",attrs:{href:"#quality-control","aria-hidden":"true"}},[e._v("#")]),e._v(" Quality Control")]),e._v(" "),n("ul",[n("li",[n("strong",[e._v("FASTQ QC metrics")]),e._v(" are generated using "),n("a",{attrs:{href:"https://github.com/OpenGene/fastp",target:"_blank",rel:"noopener noreferrer"}},[e._v("fastp"),n("OutboundLink")],1),e._v(".")]),e._v(" "),n("li",[n("strong",[e._v("BAM file QC metrics")]),e._v(" are generated using "),n("a",{attrs:{href:"https://github.com/tobiasrausch/alfred",target:"_blank",rel:"noopener noreferrer"}},[e._v("Alfred"),n("OutboundLink")],1),e._v(".")]),e._v(" "),n("li",[n("strong",[e._v("Hybridisation-selection metrics")]),e._v(" are generated using "),n("a",{attrs:{href:"https://software.broadinstitute.org/gatk/documentation/tooldocs/4.beta.6/picard_analysis_directed_CollectHsMetrics.php",target:"_blank",rel:"noopener noreferrer"}},[e._v("CollectHsMetrics"),n("OutboundLink")],1),e._v(". Only for exomes.")]),e._v(" "),n("li",[n("strong",[e._v("Contamination and concordance metrics")]),e._v(" for tumor-normal pairs using "),n("a",{attrs:{href:"https://github.com/mskcc/Conpair",target:"_blank",rel:"noopener noreferrer"}},[e._v("Conpair"),n("OutboundLink")],1),e._v(".")])])])},[],!1,null,null,null);t.default=a.exports}}]);