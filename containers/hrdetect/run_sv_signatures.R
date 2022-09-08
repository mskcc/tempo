#!/usr/bin/env Rscript

suppressPackageStartupMessages({
library(signature.tools.lib)
library(data.table)
library(tools)
library(stringr)
library(BSgenome.Hsapiens.1000genomes.hs37d5)
library(BSgenome.Hsapiens.UCSC.hg38)
library(plyr)
library(dplyr)
library(getopt)
})

how_to <- function(){
    message(" ")
    message("This script runs the signature calling function for SVs from signature.tools.lib")
    message("Run this script as follows:")
    message(" ")
    message("run_sv_signatures.R -i BEDPE -g GENOME -s NAME -o OUTDIR -n NPARALLEL")
    message(" ")
    message("    -i BEDPE	Input bedpe with candidate SV calls")
    message("    -g GENOME	hg19/hg38")
    message("    -s NAME	Sample name")
    message("    -o OUTDIR	Output directory")
    message("    -n NPARALLEL	Number of threads")
    message("    -h		Show this explanation")
}

spec = matrix(c(
		'input',	'i', 1, 'character',
		'genome',	'g', 1, 'character',
		'sampleName',	's', 1, 'character',
		'outDir',	'o', 2, 'character',
		'nparallel',	'n', 2, 'integer',
		'help',	'h', 0, 'logical'
		), byrow=TRUE, ncol=4)
opt = getopt(spec)

if ( !is.null(opt$help) ) {
	how_to()
	q(status=1,save = "no")
}

input.file <- opt$input
genome.v = opt$genome
if (genome.v == 'hg19'){
	ref.genome = BSgenome.Hsapiens.1000genomes.hs37d5
} else if (genome.v == 'hg38'){
	ref.genome = BSgenome.Hsapiens.NCBI.GRCh38
} else {
	stop("Input genome must be hg19/hg38")
}

sampleName <- opt$sampleName
outDir <- opt$outDir
nparallel <- ifelse( is.null(opt$nparallel),1,opt$nparallel)
message("## creating output directory if it does not exist ##")
if ( is.null(opt$outDir)){ outDir <- "./" } else {outDir <- opt$outDir}
if ( ! dir.exists(outDir)){ dir.create(outDir) }

message("## read and pre-process input ##")
this_sv <- fread(cmd = paste("grep -v '^##'", input.file), header = T, data.table = F)
if ("FILTER" %in% names(this_sv)){ this_sv <- this_sv[this_sv$FILTER=="PASS",] }
this_sv$sample <- sampleName
setnames(this_sv, "#CHROM_A", "chrom1")
setnames(this_sv, "START_A", "start1")
setnames(this_sv, "END_A", "end1")
setnames(this_sv, "CHROM_B", "chrom2")
setnames(this_sv, "START_B", "start2")
setnames(this_sv, "END_B", "end2")
setnames(this_sv, "STRAND_A", "strand1")
setnames(this_sv, "STRAND_B", "strand2")
svclass_dict <- c("BND"="translocation","INV"="inversion","DEL"="deletion","DUP"="tandem-duplication")
this_sv$svclass <- stringr::str_replace_all(string  = this_sv$TYPE,
					    pattern = svclass_dict)
this_sv <- this_sv[,unlist(strsplit("chrom1,start1,end1,chrom2,start2,end2,strand1,strand2,sample,svclass",split=","))]

SV_bedpe_file <- paste(outDir,paste0(sampleName,".reformat.bedpe"),sep="/")
writeTable(this_sv, SV_bedpe_file)
names(SV_bedpe_file) <- sampleName

message("## use signature.tools.lib functions ##")
message(paste0("# Version: ",as.character(packageVersion("signature.tools.lib"))))
cat_sv <- bedpeToRearrCatalogue(this_sv)
randomSeed <- NULL
set.seed(randomSeed)

sig.arglist = list()
sig.arglist[["genome.v"]] = genome.v
sig.arglist[["nparallel"]] = nparallel
if (strsplit(as.character(packageVersion("signature.tools.lib")),"\\.")[[1]][1] %in% c("0","1")){
	sig.arglist[["cat"]] = cat_sv
	sig.arglist[["signature_data_matrix"]] = signature.tools.lib:::RS.Breast560
	res <- do.call(SignatureFit_withBootstrap, sig.arglist)
	plotRearrSignatures(res$cat,  output_file = paste(outDir,paste0(sampleName,"_catalogues.pdf"),sep="/"))
	writeTable(t(res$E_median_filtered),paste(outDir,paste0(sampleName,"_exposures.tsv"),sep="/"))
} else {

	sig.arglist[["useBootstrap"]] = TRUE
	sig.arglist[["fit_method"]] = "Fit" 
	sig.arglist[["SV_bedpe_files"]] = SV_bedpe_file
	sig.arglist[["signature_version"]] = "RefSigv2"
	res <- do.call(signatureFit_pipeline, sig.arglist)
	plotSignatures(res$catalogues, output_file = paste(outDir,paste0(sampleName,"_catalogues.pdf"),sep="/"),ncolumns=1)
	writeTable(res$fitResults$exposures, paste(outDir,paste0(sampleName,"_exposures.tsv"),sep="/"))
}
saveRDS(res,file="result.rds")
 
