#!/usr/bin/env Rscript

# Prepare sample statistics for BRASS SV caller
# Outputs <tag>.samplestatistics.txt, <tag>.facets.copynumber.csv, <tag>.facets.filtered.copynumber.csv
# Usage: Rscript generate_samplestatistics.R <rdatafile> <tag> 
args <- commandArgs(trailingOnly=TRUE)
rdatafile <- args[1]
tag <- args[2]

ssfile <- paste0(tag,".samplestatistics.txt") 
cnvfile <- paste0(tag,".facets.copynumber.csv") 
cnvfile.filter <- paste0(tag,".facets.filtered.copynumber.csv") 

load(rdatafile)
slice <- out$out[out$out$chrom==23,]
if (sum(slice$nhet)/sum(slice$num.mark) < .01 ){GenderChrFound <- 'Y'} else {GenderChrFound <- 'N' }
cat(paste('Ploidy',fit$ploidy,"\n"), file=ssfile)
cat(paste('rho',fit$purity,"\n"), file=ssfile,append = T) # rho = purity
cat(paste('GenderChr','Y',"\n"), file=ssfile,append = T)
cat(paste('GenderChrFound',GenderChrFound,"\n"), file=ssfile,append = T)

cncf <- fit$cncf

cncf$n.tcn <- ifelse(
	cncf$chrom == 24,
	ifelse(GenderChrFound == 'Y', 1, 0),
	ifelse(cncf$chrom < 23, 2, ifelse(GenderChrFound == 'Y', 1, 2))
)
cncf$n.lcn <- ifelse(
	cncf$chrom == 24,
	0,
	ifelse(cncf$chrom < 23, 1, ifelse(GenderChrFound == 'Y', 0, 1))
)

cncf.reorder <- cncf[,c(unlist(strsplit("seg,chrom,start,end,n.tcn,n.lcn,tcn,lcn",",")))]

filter.na.vector <- function(x){
	any(is.na(x))
}
cncf.reorder.filter <- cncf.reorder[!apply(cncf.reorder,1,filter.na.vector),]

write.table(
	cncf.reorder, 
	file=cnvfile, 
	sep=",", 
	row.names=F, 
	col.names=F, 
	quote=F
)
write.table(
	cncf.reorder.filter, 
	file=cnvfile.filter, 
	sep=",", 
	row.names=F, 
	col.names=F, 
	quote=F
)
