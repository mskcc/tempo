#!/usr/bin/env Rscript

## Tabulate alignment metrics from Alfred and CollectHsMetrics across multiple samples
## author  = Philip Jonsson
## email   = jonssonp@mskcc.org
## version = 0.1.0
## status  = Dev

suppressPackageStartupMessages({
    library(data.table)
    library(dplyr)
    library(parallel)
})

args = commandArgs(TRUE)

if (is.null(args) | length(args)<1) {
    message('Usage: create-aggregate-qc-file.R assay\n')
    quit()
}

# Read Alfred output ----------------------------------------------------------------------------------------------

read_alfred_qc = function(file) {
    fread(cmd = paste('zgrep ^ME', file, '| cut -f 2-')) %>% 
        mutate(TotalReads = round(`#Mapped`/MappedFraction)) %>% 
        select(Sample,
               TotalReads,
               NumberDuplicateMarked = `#DuplicateMarked`,
               FractionDuplicateMarked = DuplicateFraction,
               NumberMapped = `#Mapped`,
               FractionMapped = MappedFraction,
               NumberUnmapped = `#Unmapped`,
               FractionUnmapped = UnmappedFraction,
               NumberSecondaryAlignments = `#SecondaryAlignments`,
               FractionSecondaryAlignments = SecondaryAlignmentFraction,
               TotalPairs = `#Pairs`,
               NumberMappedProperPairs = `#MappedProperPair`,
               FractionMappedProperPairs = MappedProperFraction,
               MedianReadLength,
               MedianInsertSize,
               MedianCoverage)
}

# Read CollectHsMetrics output ------------------------------------------------------------------------------------

read_hs_metrics = function(file) {
    fread(cmd = paste('grep -A2 \"## METRICS\"', file, '| tail -n +2')) %>% 
        mutate(Sample = gsub('.hs_metrics.txt', '', file)) %>% 
        select(
            Sample,
            BasesOnBait = ON_BAIT_BASES,
            BasesOnTarget = ON_TARGET_BASES,
            MeanBaitCoverage = MEAN_BAIT_COVERAGE,
            MeanTargetCoverage = MEAN_TARGET_COVERAGE,
            MedianTargetCoverage = MEDIAN_TARGET_COVERAGE,
            MaxTargetCoverage = MAX_TARGET_COVERAGE,
            FractionTargetsZeroCoverage = ZERO_CVG_TARGETS_PCT,
            FractionTargets10X = PCT_TARGET_BASES_10X,
            FractionTargets20X = PCT_TARGET_BASES_20X,
            FractionTargets50X = PCT_TARGET_BASES_50X,
            FractionTargets100X = PCT_TARGET_BASES_100X)
}

# Execute in run directory ----------------------------------------------------------------------------------------

if (!interactive()) {
    
    # Run-time parameters
    num_cores = detectCores()
    assay_type = args[1]
    
    # Parse output from Alfred 
    alfred_files = dir(getwd(), pattern = '*.alfred.tsv.gz$')
    outMatrix = mclapply(alfred_files, read_alfred_qc, mc.cores = num_cores)
    outMatrix = do.call('rbind', outMatrix)
    
    # For exomes, parse output from CollectHsMetrics
    if (assay_type == 'wes') {
        hs_metrics_files = dir(getwd(), pattern = '*.hs_metrics.txt$')
        hs_metrics_out = mclapply(hs_metrics_files, read_hs_metrics, mc.cores = num_cores) 
        hs_metrics_out = do.call('rbind', hs_metrics_out)
        
        outMatrix = left_join(outMatrix, hs_metrics_out, by = 'Sample') %>% 
            select(-MedianCoverage)
    }
    
    fwrite(outMatrix, 'alignment_qc.tsv', sep = '\t')

}
