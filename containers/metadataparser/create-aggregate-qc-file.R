#!/usr/bin/env Rscript

## Tabulate alignment metrics from Alfred and CollectHsMetrics across multiple samples
## author  = Philip Jonsson
## contributor = Evan Biederstedt, evan.biederstedt@gmail.com, 2019 ; Anne Marie Noronha, noronhaa@mskcc.org, 2020

## email   = jonssonp@mskcc.org, evan.biederstedt@gmail.com, noronhaa@mskcc.org
## version = 0.2.1
## status  = Dev

suppressPackageStartupMessages({
    library(data.table)
    library(parallel)
})

args = commandArgs(TRUE)

if (is.null(args) | length(args)<1) {
    message('Usage: create-aggregate-qc-file.R assay\n')
    quit()
}

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


## Execute in run directory ----------------------------------------------------------------------------------------

if (!interactive()) {
    
    # Run-time parameters
    assay_type = args[1]
    num_cores = args[2]
    
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
    }
    
    fwrite(outMatrix, 'alignment_qc.txt', sep = '\t')

}

