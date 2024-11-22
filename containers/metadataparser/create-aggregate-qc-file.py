#!/usr/bin/env python
# -*- coding: utf-8 -*-
 
"""
Custom aggregation of data from Alfred and GATK CollectHSMetrics.
Usage: filter-vcf.py assay
Output: 'input_filename.filter.vcf'
"""

__author__  = "Philip Jonsson"
__email__   = "jonssonp@mskcc.org"
__version__ = "0.1.0"
__status__  = "Dev"

import sys, glob, subprocess, multiprocessing
import pandas as pd
from io import StringIO
from joblib import Parallel, delayed

numCores = multiprocessing.cpu_count()
assayType = sys.argv[1]

def main():
    alfredFiles = glob.glob('*.alfred.tsv.gz')

    alfredOut = Parallel(n_jobs = numCores)(delayed(readAlfredQC)(fl) for fl in alfredFiles)
    outMatrix = pd.concat(alfredOut)

    if assayType == "wes":
        hsMetricsFiles = glob.glob('*.hs_metrics.txt')
        hsMetricsOut = Parallel(n_jobs = numCores)(delayed(readHsMetrics)(fl) for fl in hsMetricsFiles)
    
        hsMetricsOutMatrix = pd.concat(hsMetricsOut)
        outMatrix = pd.merge(outMatrix, hsMetricsOutMatrix, on = "Sample")

        outMatrix = outMatrix.drop(columns = "MedianCoverage")

    outMatrix.to_csv("alignment_qc.tsv", sep = "\t", index = False)

def readAlfredQC(file):
    readCmd = "zgrep ^ME %s | cut -f 2-" % file
    proc = subprocess.Popen(readCmd, shell = True, stdout = subprocess.PIPE)
    a = StringIO(proc.communicate()[0].decode('utf-8'))
    b = pd.read_csv(a, sep = "\t")
    b["TotalReads"] = round(b["#Mapped"]/b["MappedFraction"])

    c = b.rename(columns = {
        "#DuplicateMarked": "NumberDuplicateMarked",
        "DuplicateFraction": "FractionDuplicateMarked",
        "#Mapped": "NumberMapped",
        "MappedFraction": "FractionMapped",
        "#Unmapped": "NumberUnmapped",
        "UnmappedFraction": "FractionUnmapped",
        "#SecondaryAlignments": "NumberSecondaryAlignments",
        "SecondaryAlignmentFraction": "FractionSecondaryAlignments",
        "#Pairs": "TotalPairs",
        "#MappedProperPair": "NumberMappedProperPairs",
        "MappedProperFraction": "FractionMappedProperPairs"
    })
        
    return c[["Sample", "TotalReads", "NumberDuplicateMarked", "FractionDuplicateMarked", "NumberMapped", "FractionMapped", "NumberUnmapped", "FractionUnmapped", "NumberSecondaryAlignments", "FractionSecondaryAlignments", "TotalPairs", "NumberMappedProperPairs", "FractionMappedProperPairs", "MedianReadLength", "MedianInsertSize", "MedianCoverage"]]

def readHsMetrics(file):
    readCmd = "grep -A2 \"## METRICS\" %s | tail -n +2" % file
    proc = subprocess.Popen(readCmd, shell = True, stdout = subprocess.PIPE)
    a = StringIO(proc.communicate()[0].decode('utf-8'))
    b = pd.read_csv(a, sep = "\t")

    b['Sample'] = file.split(".hs_metrics.txt")[0]

    c = b.rename(columns = {
        "ON_BAIT_BASES": "BasesOnBait",
        "ON_TARGET_BASES": "BasesOnTarget",
        "MEAN_BAIT_COVERAGE": "MeanBaitCoverage",
        "MEAN_TARGET_COVERAGE": "MeanTargetCoverage",
        "MEDIAN_TARGET_COVERAGE": "MedianTargetCoverage",
        "MAX_TARGET_COVERAGE": "MaxTargetCoverage",
        "ZERO_CVG_TARGETS_PCT": "FractionTargetsZeroCoverage",
        "PCT_TARGET_BASES_10X": "FractionTargets10X",
        "PCT_TARGET_BASES_20X": "FractionTargets20X",
        "PCT_TARGET_BASES_50X": "FractionTargets50X",
        "PCT_TARGET_BASES_100X": "FractionTargets100X"
    })

    return c[["Sample", "BasesOnBait", "BasesOnTarget", "MeanBaitCoverage", "MeanTargetCoverage", "MedianTargetCoverage", "MaxTargetCoverage", "FractionTargetsZeroCoverage", "FractionTargets10X", "FractionTargets20X", "FractionTargets50X", "FractionTargets100X"]]

if __name__ == "__main__":
    main()
