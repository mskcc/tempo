#!/usr/bin/env Rscript

# __author__      = "Sam Tischfield"
# __email__       = "tischfis@mskcc.org"
# __contributor__ = "Anne Marie Noronha (noronhaa@mskcc.org)"
# __version__     = "0.0.1"
# __status__      = "Dev"

suppressPackageStartupMessages({
	library(BioCircos)
	library(data.table)
	library(tidyverse)
	library(getopt)
})

how_to <- function(){
    message(" ")
    message("This script runs the signature calling function for SVs from signature.tools.lib")
    message("Run this script as follows:")
    message(" ")
    message("generate_biocircos.R -i BEDPE -g GENOME -s NAME -o OUTDIR -n NPARALLEL")
    message(" ")
    message("    -b BEDPE	Input bedpe with candidate SV calls")
	message("    -c CNCF	Input CNCF table from Facets")
    message("    -g GENOME	hg19/hg38 (default hg19)")
    message("    -s NAME	Sample name")
    message("    -h		Show this explanation")
}

spec = matrix(c(
		'bedpe.path',	'b', 1, 'character',
		'cncf.path',	'c', 1, 'character',
		'genome',	'g', 2, 'character',
		'sampleName',	's', 1, 'character',
		'help',	'h', 0, 'logical'
		), byrow=TRUE, ncol=4)
opt = getopt(spec)
if ( !is.null(opt$help) ) {
	how_to()
	q(status=1,save = "no")
}

if ( !is.null(opt$genome) ) {
	genome_v <- 'hg19'
} else if (opt$genome %in% c('hg19','hg38')){
	genome_v <- opt$genome
} else {
	how_to()
	q(status=1,save = "no")
}

length.connection = pipe(paste("cat ", opt$bedpe.path, " | grep \"^##\" | wc -l", sep = ""))
header.length = as.numeric(trimws(readLines(con = length.connection, n = 1)))
close(length.connection)
bedpe=fread(opt$bedpe.path,
            skip=header.length)
cnv=fread(opt$cncf.path)
if (!"tcn.em" %in% names(cnv)){
	# assume ascat format
	names(cnv) <- unlist(str_split("index,chrom,loc.start,loc.end,tcn.em,lcn.em,normal_tcn,normal_lcn", pattern=","))
} 
cnv %>%
	select(c(chrom,loc.start,loc.end,tcn.em,lcn.em)) %>%
	dplyr::rename(tcn = 'tcn.em') %>%
	dplyr::rename(lcn = 'lcn.em') -> cnv



plot.biocircos <- function(bedpe,cnv,sample_name,genome_v){
  # initiate variable to store tracks
  tracklist <- list()

  # add bedpe data
  bedpe %>%
    dplyr::rename(CHROM_A=`#CHROM_A`) -> bedpe
  gen.BioCircosLinkTrack <- function(trackname,bedpe,color){
    b <- BioCircosLinkTrack(trackname,
                            #labels=paste(trackname, bedpe$START_A, bedpe$END_A, sep=":"),
                            labels=bedpe$ID,
                            gene1Starts = bedpe$START_A,
                            gene1Ends = bedpe$END_A,
                            gene1Chromosomes = bedpe$CHROM_A,
                            gene2Chromosomes = bedpe$CHROM_B,
                            gene2Starts=bedpe$START_B,
                            gene2Ends = bedpe$END_B,
                            maxRadius = .5,
                            color = color,
                            displayLabel = F
    )
    }
  for (i in list(c("BND","black"),
                 c("DEL","red"),
                 c("DUP","green"),
                 c("INV","blue"),
                 c("INS","yellow")
                 )
       ){
    bedpe.slice <- bedpe %>% filter(TYPE==i[1])
    if (dim(bedpe.slice)[[1]] > 0){
    tracklist[[i[1]]] = gen.BioCircosLinkTrack(i[1], bedpe.slice, i[2])
    }
  }

  # add cnv data
  # convert 23 to X and 24 to Y
  print(head(cnv))
  cnv %>%
    mutate(chrom=gsub(pattern = 23,replacement = "X",x = chrom)) %>%
    mutate(chrom=gsub(pattern = 24,replacement = "Y",x = chrom)) %>%
    mutate(tcn.em=ifelse(tcn>=10,10,tcn)) ->cnv_test
  cnv.range=c(0,max(cnv_test$tcn.em,na.rm = T))
  cnv.range
  tracklist$tcn.em =
    BioCircosCNVTrack('cnv_track',
                      chromosomes = cnv_test$chrom,
                      starts = cnv_test$loc.start,
                      ends = cnv_test$loc.end,
                      values = cnv_test$tcn,
                      color = "black",
                      range = cnv.range)
  tracklist$lcn.em =
    BioCircosCNVTrack('cnv_track',
                      chromosomes = cnv_test$chrom,
                      starts = cnv_test$loc.start,
                      ends = cnv_test$loc.end,
                      values = cnv_test$lcn,
                      color = "red",
                      range = cnv.range)

  tracklist$background = BioCircosBackgroundTrack("arcs_background", colors = "#2222EE")
  
  
  tracklist <- Reduce('+', tracklist)
  plot=BioCircos(tracklist,
                 genomeFillColor = "PuOr",
                 genome = genome_v,
                 chrPad = 0.02,
                 displayGenomeBorder = T,
                 yChr =  FALSE,
                 genomeTicksDisplay = FALSE,
                 genomeLabelTextSize = "8pt",
                 genomeLabelDy = 0)
  return(plot)
}

plot <- plot.biocircos(bedpe=bedpe,cnv=cnv,sample_name = opt$sampleName, genome_v=genome_v )
htmlwidgets::saveWidget(plot, paste0(opt$sampleName,".circos.html"), selfcontained = T, libdir = "lib")
