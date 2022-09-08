#!/usr/bin/env Rscript


# __author__      = "Vasilisa Rudneva"
# __email__       = "rudnevav@mskcc.org"
# __contributor__ = "Yixiao Gong (gongy@mskcc.org); Anne Marie Noronha (noronhaa@mskcc.org)"
# __version__     = "1.0"


# ## Installation:
# install.packages("devtools")
# BiocManager::install("BSgenome")
# BiocManager::install("RcppProgress")
# BiocManager::install("BSgenome.Mmusculus.UCSC.mm10")
# BiocManager::install("BSgenome.Hsapiens.1000genomes.hs37d5")
# BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
# install.packages("https://cran.r-project.org/src/contrib/Archive/NNLM/NNLM_0.4.3.tar.gz", repo=NULL)
# #devtools::install_github('linxihui/NNLM')
# devtools::install_github('Nik-Zainal-Group/signature.tools.lib')

suppressPackageStartupMessages({
library(signature.tools.lib)
library(data.table)
library(tools)
library(stringr)
library(BSgenome.Hsapiens.1000genomes.hs37d5)
library(BSgenome.Hsapiens.UCSC.hg38)
library(plyr)
library(dplyr)
})

args = commandArgs(trailingOnly=TRUE)


inputTSV = args[1]
genome_version = args[2]
if (genome_version == 'hg19'){
	ref.genome = BSgenome.Hsapiens.1000genomes.hs37d5
} else if (genome_version == 'hg38'){
	ref.genome = BSgenome.Hsapiens.NCBI.GRCh38
} else {
	stop("Input genome must be hg19/hg38")
}

n_parallel = as.integer(args[3])
if(is.na(n_parallel)){n_parallel=1}

message(paste("========= Genome Version ", genome_version, "=========", sep=" "))
input.files<-fread(inputTSV, header = T, data.table = T)

sample_names<-input.files$sample

dir.create("tmp")

# 1. SV_bedpe_files
#The files should contain a header in the first line with the following columns: "chrom1", "start1", "end1", "chrom2", "start2", "end2" and "sample" (sample name). 
#In addition, either two columns indicating the strands of the mates, "strand1" (+ or -) and "strand2" (+ or -), 
#or one column indicating the structural variant class, "svclass": translocation, inversion, deletion, tandem-duplication. 
#The column "svclass" should correspond to (Sanger BRASS convention): inversion (strands +/- or -/+ and mates on the same chromosome), 
#deletion (strands +/+ and mates on the same chromosome), tandem-duplication (strands -/- and mates on the same chromosome), 
#translocation (mates are on different chromosomes)..

correctSV <- function(this_sample){
  cat(this_sample)
  this_sv<-fread(cmd = paste("grep -v '^##'", input.files[input.files$sample==this_sample,]$sv), header = T, data.table = F)
  this_sv<-this_sv[this_sv$FILTER=="PASS",]
  setnames(this_sv, "#CHROM_A", "chrom1")
  setnames(this_sv, "START_A", "start1")
  setnames(this_sv, "END_A", "end1")
  setnames(this_sv, "CHROM_B", "chrom2")
  setnames(this_sv, "START_B", "start2")
  setnames(this_sv, "END_B", "end2")
  setnames(this_sv, "STRAND_A", "strand1")
  setnames(this_sv, "STRAND_B", "strand2")
  svclass_dict <- c("BND"="translocation","INV"="inversion","DEL"="deletion","DUP"="tandem-duplication") 
  this_sv$svclass <- stringr::str_replace_all(string = this_sv$TYPE,
					      pattern= svclass_dict)
  print(table(this_sv$svclass))
  this_sv$sample=this_sample
  this_sv <- this_sv[,unlist(strsplit("chrom1,start1,end1,chrom2,start2,end2,strand1,strand2,sample,svclass",split=","))]
  filename<-paste0("tmp/",this_sample, ".sv")
  write.table(file = filename, x = this_sv, quote = F, row.names = F, col.names = T, sep = "\t")
}
message("")
message("========= Preprocessing SV inputs =========")
invisible(lapply(sample_names,correctSV))
SV_bedpe_files <- paste0("tmp/",sample_names,".sv")
names(SV_bedpe_files) <- sample_names


# 2 Indels_tab_files and SNV_tab_files
#list of file names corresponding to Indels/SNV TAB files to be used to classify Indels and compute the proportion of indels at micro-homology or 96-channel substitution catalogues. 
#This should be a named vector, where the names indicate the sample name, 
#so that each file can be matched to the corresponding row in the data_matrix input. 
#The files should only contain indels (or SNV) and should already be filtered according to the user preference, 
#as all indels in the file will be used and no filter will be applied. 
#Each File contains indels from a single sample and the following minimal columns: chr, position, REF, ALT.
correctMutations <- function(this_sample){
  cat(this_sample)
  this_mutations<-fread(input.files[input.files$sample==this_sample,]$mutations, header = T, data.table = F)
  setnames(this_mutations, "Chromosome", "chr")
  setnames(this_mutations, "vcf_pos", "position")
  setnames(this_mutations, "vcf_id", "ID")
  setnames(this_mutations, "Reference_Allele", "REF")
  setnames(this_mutations, "Allele", "ALT")
  setnames(this_mutations, "vcf_qual", "QUAL")
  this_indels<-this_mutations[this_mutations$Variant_Type %in% c("DEL", "DNP", "INS", "TNP"),]
  this_indels <- this_indels %>%
      dplyr::rowwise() %>%
      mutate(left_b = as.character(ref.genome[[chr]][position])) %>%
      mutate(REF = ifelse(ALT == "-",paste0(left_b,REF),REF),
	     ALT = ifelse(ALT == "-",left_b,ALT)) %>%
      mutate(ALT = ifelse(REF == "-",paste0(left_b,ALT),ALT),
	     REF = ifelse(REF == "-",left_b,REF)) %>%
      select(-c(left_b))

  this_snv<-this_mutations[!this_mutations$Variant_Type %in% c("DEL", "DNP", "INS", "TNP"),]
  print(table(rbind(this_indels, this_snv)$Variant_Type))
  write.table(file = paste0("tmp/", this_sample, ".indels"), x = this_indels, quote = F, row.names = F, col.names = T, sep = "\t")
  write.table(file = paste0("tmp/", this_sample, ".snv"), x = this_snv, quote = F, row.names = F, col.names = T, sep = "\t")
}
message("")
message("========= Preprocessing Indels and SNV inputs =========")
invisible(lapply(sample_names,correctMutations))
Indels_tab_files <- paste0("tmp/",sample_names,".indels")
SNV_tab_files <- paste0("tmp/",sample_names,".snv")
names(Indels_tab_files) <- sample_names
names(SNV_tab_files) <- sample_names

# 3. CNV_tab_files
# list of file names corresponding to CNV TAB files (similar to ASCAT format) 
#to be used to compute the HRD-LOH index. This should be a named vector, 
#where the names indicate the sample name, 
#so that each file can be matched to the corresponding row in the data_matrix input. 
#The files should contain a header in the first line with the following columns: 
#'seg_no', 'Chromosome', 'chromStart', 'chromEnd', 
#'total.copy.number.inNormal', 'minor.copy.number.inNormal', 'total.copy.number.inTumour', 'minor.copy.number.inTumour'
correctCNV <- function(this_sample){
  cat(this_sample)
  this_cnv<-fread(input.files[input.files$sample==this_sample,]$cnv, header = F, data.table = F)
  # Input file does not have header. Be alert. https://github.com/cancerit/ascatNgs/wiki/Protocol-Correction---format-of-copynumber.caveman.csv
  colnames(this_cnv)<-c("seg_no", "Chromosome", "chromStart", "chromEnd", "total.copy.number.inNormal", "minor.copy.number.inNormal", "total.copy.number.inTumour", "minor.copy.number.inTumour")
  print(table(this_cnv$Chromosome))
  filename<-paste0("tmp/", this_sample, ".cnv")
  write.table(file = filename, x = this_cnv, quote = F, row.names = F, col.names = T, sep = "\t")
}
message("")
message("========= Preprocessing CNV inputs =========")
invisible(lapply(sample_names,correctCNV))
CNV_tab_files <- paste0("tmp/",sample_names,".cnv")
names(CNV_tab_files) <- sample_names


#load SNV data and convert to SNV mutational catalogues
SNVcat_list <- list()
message("")
message("========= Converting SNV mutational catalogues =========")
for (i in 1:length(SNV_tab_files)){
  message(sample_names[i])
  tmpSNVtab <- read.table(SNV_tab_files[i],sep = "\t", fill = T, quote="", header = TRUE,check.names = FALSE, stringsAsFactors = FALSE)
  res <- tabToSNVcatalogue(subs = tmpSNVtab,genome.v = genome_version)
  colnames(res$catalogue) <- sample_names[i]
  SNVcat_list[[i]] <- res$catalogue
}
SNV_catalogues <- do.call(cbind,SNVcat_list)

#Initialize feature matrix
col_hrdetect <- c("del.mh.prop", "SNV3", "SV3", "SV5", "hrd", "SNV8")
input_matrix <- matrix(NA,nrow = length(sample_names), ncol = length(col_hrdetect), dimnames = list(sample_names,col_hrdetect))

# Compute the proportion of indels at micro-homology 
Indel.del.mh.prop_list <- list()
message("")
message("========= Computing the proportion of indels at micro-homology =========")
for (i in 1:length(Indels_tab_files)){
  message(sample_names[i])
  tmpIndeltab <- read.table(Indels_tab_files[i],sep = "\t", fill = T, quote="", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
  res<-tabToIndelsClassification(tmpIndeltab,sample_names[i], genome_version)
  Indel.del.mh.prop_list[[i]] <- res$count_proportion
}
names(Indel.del.mh.prop_list)<-sample_names
Indel.del.mh.prop <- do.call(rbind,Indel.del.mh.prop_list)

input_matrix[rownames(Indel.del.mh.prop),"del.mh.prop"] <- Indel.del.mh.prop[,"del.mh.prop"]

# Run the pipeline
message("")
message("========= Running HRDetect Pipeline =========")
hrd.arglist = list()
hrd.arglist[["data_matrix"]] = input_matrix
hrd.arglist[["genome.v"]] = genome_version
hrd.arglist[["SV_bedpe_files"]] = SV_bedpe_files
hrd.arglist[["Indels_tab_files"]] = Indels_tab_files
hrd.arglist[["CNV_tab_files"]] = CNV_tab_files
hrd.arglist[["SNV_catalogues"]] = SNV_catalogues
hrd.arglist[["nparallel"]] = n_parallel
if (strsplit(as.character(packageVersion("signature.tools.lib")),"\\.")[[1]][1] %in% c("0","1")){
	hrd.arglist[["signature_type"]] = "COSMIC"
} else {
	hrd.arglist[["SNV_signature_version"]] = "COSMICv2"
}
res <- do.call(HRDetect_pipeline, hrd.arglist)
if (1==0) {
res<-HRDetect_pipeline(input_matrix,
                       genome.v = genome_version,
                       SV_bedpe_files = SV_bedpe_files, 
                       Indels_tab_files  = Indels_tab_files, 
                       CNV_tab_files = CNV_tab_files,
                       #SNV_tab_files = SNV_tab_files, 
                       SNV_catalogues = SNV_catalogues,
                       nparallel = n_parallel)
}

#save HRDetect scores
hrdetect_output = as.data.table(res$hrdetect_output, keep.rownames="sample")
write.table(hrdetect_output,file = paste0(basename(file_path_sans_ext(inputTSV)),".hrdetect.tsv"), row.names=F, quote=F, sep = "\t")

# Cleanup
unlink("tmp", recursive = T)
