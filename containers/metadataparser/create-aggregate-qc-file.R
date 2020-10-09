#!/usr/bin/env Rscript

## Tabulate alignment metrics from Alfred and CollectHsMetrics across multiple samples
## author  = Philip Jonsson
## contributor = Evan Biederstedt, evan.biederstedt@gmail.com, 2019 ; Anne Marie Noronha, noronhaa@mskcc.org, 2020

## email   = jonssonp@mskcc.org, evan.biederstedt@gmail.com, noronhaa@mskcc.org
## version = 0.2.2
## status  = Dev

suppressPackageStartupMessages({
    library(data.table)
    library(parallel)
    library(argparse)
})

#args = commandArgs(TRUE)

#if (is.null(args) | length(args)<1) {
#    message('Usage: create-aggregate-qc-file.R assay\n')
#    quit()
#} 
## default values now exist for each argument. error message is no longer necessary.

parser = ArgumentParser(description = 'Aggregate QC metrics across samples')
parser$add_argument('-n', '--num-cores', type = 'integer', required = TRUE, default = 1,
                    help = 'number of cores')
parser$add_argument('-a', '--assay-type', type = 'character', required = TRUE, default = 'wes',
                    help = 'Assay type')
args = parser$parse_args()

## Parse Alfred output ----------------------------------------------------------------------------------------------
read_alfred_qc = function(file) {
    input_file = fread(cmd = paste('zgrep ^ME', file, '| cut -f 2-')) 
    ## define column TotalReads
    input_file[, TotalReads := round(`#Mapped`/MappedFraction)]
    ##
    ## only select columns
    ## c('Sample', 'TotalReads','#DuplicateMarked','DuplicateFraction','#Mapped','MappedFraction','#Unmapped',
    ##     'UnmappedFraction','#SecondaryAlignments','SecondaryAlignmentFraction','#Pairs','#MappedProperPair','MappedProperFraction','MedianReadLength','MedianInsertSize','MedianCoverage')
    ## 
    ## For clarity, doing this in two steps
    ## subset columns
    input_file = input_file[, c('Sample', 'TotalReads','#DuplicateMarked','DuplicateFraction','#Mapped','MappedFraction',
          '#Unmapped','UnmappedFraction','#SecondaryAlignments','SecondaryAlignmentFraction','#Pairs',
          '#MappedProperPair','MappedProperFraction','MedianReadLength','MedianInsertSize','MedianCoverage')]

    setnames(input_file, 
        old = c('#DuplicateMarked', 'DuplicateFraction', '#Mapped', 'MappedFraction', '#Unmapped', 'UnmappedFraction', '#SecondaryAlignments', 'SecondaryAlignmentFraction', '#Pairs','#MappedProperPair', 'MappedProperFraction'), 
        new = c('NumberDuplicateMarked', 'FractionDuplicateMarked', 'NumberMapped', 'FractionMapped', 'NumberUnmapped', 'FractionUnmapped', 'NumberSecondaryAlignments', 'FractionSecondaryAlignments', 'TotalPairs', 'NumberMappedProperPairs', 'FractionMappedProperPairs'))
   
    return(input_file)
} 


## Read CollectHsMetrics output ------------------------------------------------------------------------------------
read_hs_metrics = function(file) {
    input_file = fread(cmd = paste('grep -A2 \"## METRICS\"', file, '| tail -n +2')) 
    ### define column Sample
    input_file[, Sample := gsub('.hs_metrics.txt', '', file)]
    ## For clarity, doing this in two steps
    ## subset columns
    input_file = input_file[, c("Sample", "ON_BAIT_BASES", "ON_TARGET_BASES", "MEAN_BAIT_COVERAGE", "MEAN_TARGET_COVERAGE", "MEDIAN_TARGET_COVERAGE", "MAX_TARGET_COVERAGE",
        "ZERO_CVG_TARGETS_PCT", "PCT_TARGET_BASES_10X", "PCT_TARGET_BASES_20X", "PCT_TARGET_BASES_50X", "PCT_TARGET_BASES_100X")]

    setnames(input_file, 
        old = c("Sample", "ON_BAIT_BASES", "ON_TARGET_BASES", "MEAN_BAIT_COVERAGE", "MEAN_TARGET_COVERAGE", "MEDIAN_TARGET_COVERAGE", "MAX_TARGET_COVERAGE", "ZERO_CVG_TARGETS_PCT", "PCT_TARGET_BASES_10X", "PCT_TARGET_BASES_20X", "PCT_TARGET_BASES_50X", "PCT_TARGET_BASES_100X"), 
        new = c("Sample", "BasesOnBait", "BasesOnTarget", "MeanBaitCoverage","MeanTargetCoverage","MedianTargetCoverage", "MaxTargetCoverage", "FractionTargetsZeroCoverage", "FractionTargets10X", "FractionTargets20X", "FractionTargets50X", "FractionTargets100X"))

    return(input_file)
}

read_wgs_metrics = function(file) {
    input_file = fread(cmd = paste('grep -A2 \"## METRICS\"', file, '| tail -n +2')) 
    input_file[, Sample := gsub('.wgs_metrics.txt', '', file)]
    input_file = input_file[,c("Sample","MEAN_COVERAGE","SD_COVERAGE","MEDIAN_COVERAGE","PCT_EXC_TOTAL","PCT_1X","PCT_10X","PCT_20X","PCT_50X","PCT_100X")]
    
    setnames(input_file, 
        old = c("Sample","MEAN_COVERAGE","SD_COVERAGE",               "MEDIAN_COVERAGE","PCT_EXC_TOTAL",        "PCT_1X","PCT_10X","PCT_20X","PCT_50X","PCT_100X"), 
        new = c("Sample","MeanCoverage", "StandardDeviationOfCoverage","MedianCoverage","FractionExcludedBases","FractionBases1X","FractionBases10X","FractionBases20X","FractionBases50X","FractionBases100X"))

    return(input_file)
}

## Execute in run directory ----------------------------------------------------------------------------------------

if (!interactive()) {
    
    # Run-time parameters
    num_cores = args$num_cores
    assay_type = args$assay_type
    
    # Parse output from Alfred 
    alfred_files = dir(getwd(), pattern = '*.alfred.tsv.gz$')
    outMatrix = mclapply(alfred_files, read_alfred_qc, mc.cores = num_cores)
    outMatrix = as.data.table(rbindlist(outMatrix, fill=TRUE)) ## same as do.call('rbind', ... )
    
    # For exomes, parse output from CollectHsMetrics
    if (assay_type == 'wes') {
        hs_metrics_files = dir(getwd(), pattern = '*.hs_metrics.txt$')
        hs_metrics_out = mclapply(hs_metrics_files, read_hs_metrics, mc.cores = num_cores) 
        hs_metrics_out = as.data.table(rbindlist(hs_metrics_out, fill=TRUE))   ## same as do.call('rbind', ... )
        setkey(outMatrix, Sample)
        setkey(hs_metrics_out, Sample)
        outMatrix = as.data.table(merge(outMatrix, hs_metrics_out, by = 'Sample'))
        outMatrix[, MedianCoverage:=NULL]
    } else if (assay_type == 'wgs'){
        wgs_metrics_files = dir(getwd(), pattern = '*.wgs_metrics.txt$')
        wgs_metrics_out = mclapply(wgs_metrics_files, read_wgs_metrics, mc.cores = num_cores) 
        wgs_metrics_out = as.data.table(rbindlist(wgs_metrics_out, fill=TRUE))   ## same as do.call('rbind', ... )
        setkey(outMatrix, Sample)
        setkey(hs_metrics_out, Sample)
        outMatrix = as.data.table(merge(outMatrix, wgs_metrics_out, by = 'Sample'))
        outMatrix[, MedianCoverage:=NULL]
    }
    
    fwrite(outMatrix, 'alignment_qc.txt', sep = '\t')

}

