#!/usr/bin/env python3


__author__  = "Evan Biederstedt"
__contributor__  = "Yixiao Gong, Anne Marie Noronha"
__email__   = "biederse@mskcc.org; evan.biederstedt@gmail.com; gongy@mskcc.org; noronhaa@mskcc.org"
__version__ = "0.2.3"
__status__  = "Dev"


"""
Parse pipeline outputs into a single *tsv file, a tab-delimited *tsv of metadata

Usage: python3 create_metadata_file.py [-h] --sampleID SAMPLEID
                               --tumorID TUMORID
                               --normalID NORMALID
                               [--facetsPurity_out FACETSPURITY_OUT]
                               [--facetsQC FACETS_QC]
                               [--MSIsensor_output MSISENSOR_OUTPUT]
                               [--mutational_signatures_output MUTATIONAL_SIGNATURES_OUTPUT]
                               [--polysolver_output POLYSOLVER_OUTPUT]
                               [--MAF_input MAF_INPUT]
                               [--coding_baits_BED CODING_BAITS_BED]

    optional arguments:
      -h, --help            show this help message and exit
      --sampleID SAMPLEID   sample ID from channel, should always exist
      --facetsPurity_out FACETSPURITY_OUT
                            FACETS purity output, *_purity.out; used to measure
                            ploidy and purity
     --facetsQC FACETS_QC
                            FACETS QC output, *.qc.txt; used to measure WGD
     --MSIsensor_output MSISENSOR_OUTPUT
                            MSIsensor output, *.msisensor.tsv
     --mutational_signatures_output MUTATIONAL_SIGNATURES_OUTPUT
                            output *txt of mutational signatures
     --polysolver_output POLYSOLVER_OUTPUT
                            output *txt of HLA genotypes, winners.hla.txt
      --MAF_input MAF_INPUT
                            annotated MAF, *maf
      --coding_baits_BED CODING_BAITS_BED
                           BED of assay-specific baits intersecting coding
                           regions, ensGene.all_CODING_exons.reference.bed  

Output: '{sampleID}_metadata.tsv'

NOTE: This is meant to be run only within the WES/WGS pipeline at https://github.com/mskcc/vaporware
The script assumes no special output subdirectory, and assumes only three assay types exist: AgilentExon_51MB_b37_v3, IDT_Exome_v1_FP_b37, WGS_b37
"""

import argparse
import pandas as pd
import pybedtools
import os

parser = argparse.ArgumentParser()
parser.add_argument('--sampleID', help = 'sample ID from channel, should always exist', required = True)
parser.add_argument('--tumorID', help = 'tumor ID from channel, should always exist', required = True)
parser.add_argument('--normalID', help = 'normal ID from channel, should always exist', required = True)
parser.add_argument('--facetsPurity_out', help = 'FACETS purity output, *_purity.out; used to measure ploidy and purity', required = False)
parser.add_argument('--facetsQC', help = 'FACETS default QC output, *.qc.txt; used to measure WGD', required = False)
parser.add_argument('--MSIsensor_output', help = 'MSIsensor output, *.msisensor.tsv', required = False)
parser.add_argument('--mutational_signatures_output', help = 'output *txt of mutational signatures', required = False)
parser.add_argument('--polysolver_output', help = 'output *txt of HLA genotypes, winners.hla.txt', required = False)
parser.add_argument('--MAF_input', help = 'annotated MAF, *maf', required = False)
parser.add_argument('--coding_baits_BED', help = 'BED of assay-specific baits intersecting coding regions, ensGene.all_CODING_exons.reference.bed', required = False)


## NOTE: 
### Only three asssys/BED files are acceptable for coding_baits_BED:
### --- Agilent: AgilentExon_51MB_b37_v3_baits.coding.sorted.merged.bed
### --- IDT: IDT_Exome_v1_FP_b37_baits.coding.sorted.merged.bed
### --- WGS: b37_wgs_calling_regions.v1.coding.sorted.merged.bed

args = parser.parse_args()

sampleID = args.sampleID 
tumorID = args.tumorID
normalID = args.normalID
facetsPurityPloidy = args.facetsPurity_out
facetsQC = args.facetsQC
MSIoutput = args.MSIsensor_output
mutationalSignatures = args.mutational_signatures_output
HLAoutput = args.polysolver_output
MAF_input = args.MAF_input
coding_baits_BED = args.coding_baits_BED

## Current plan: put this in a pandas DataFrame and save as a *tsv file

results = pd.DataFrame()

## create Tumor Normal column
results = results.assign(sample=[sampleID])


if facetsPurityPloidy is not None:
    ## parse purity and ploidy from FACETS *_purity.out
    purity_out = pd.read_csv(facetsPurityPloidy, sep="\t")
    ## create purity and ploidy columns
    results["purity"] = [purity_out.iloc[15].str.split(' = ')[0][1]]
    results["ploidy"] = [purity_out.iloc[16].str.split(' = ')[0][1]]


if facetsQC is not None:
    ## parse WGD from facets-suite *.armlevel.tsv
    qc = pd.read_csv(facetsQC, sep="\t")
    ## create WGD column
    if sum(qc.wgd) > 0:
    	results["WGD_status"] = [True]
    elif sum(qc.wgd == False) > 0:
    	results["WGD_status"] = [False]
    else:
    	results["WGD_status"] = ['NA']  ## hard-coding string NA, not NaN or None

if MSIoutput is not None:
    ## parse MSIsensor output
    MSI = pd.read_csv(MSIoutput, sep="\t")
    ## create MSI columns
    results["MSI_Total_Sites"] = MSI.Total_Number_of_Sites
    results["MSI_Somatic_Sites"] = MSI.Number_of_Somatic_Sites
    results["MSIscore"] = MSI['%']

if mutationalSignatures is not None:
    ## parse mutational signatures output
    mutsig = pd.read_csv(mutationalSignatures, sep="\t")
    ## create mutational signatures columns
    results["SBS1.observed"] = mutsig['SBS1.observed']
    results["SBS1.pvalue"] = mutsig['SBS1.pvalue']
    results["SBS2.observed"] = mutsig['SBS2.observed']
    results["SBS2.pvalue"] = mutsig['SBS2.pvalue']
    results["SBS3.observed"] = mutsig['SBS3.observed']
    results["SBS3.pvalue"] = mutsig['SBS3.pvalue']
    results["SBS4.observed"] = mutsig['SBS4.observed']
    results["SBS4.pvalue"] = mutsig['SBS4.pvalue']
    results["SBS5.observed"] = mutsig['SBS5.observed']
    results["SBS5.pvalue"] = mutsig['SBS5.pvalue']
    results["SBS6.observed"] = mutsig['SBS6.observed']
    results["SBS6.pvalue"] = mutsig['SBS6.pvalue']
    results["SBS7a.observed"] = mutsig['SBS7a.observed']
    results["SBS7a.pvalue"] = mutsig['SBS7a.pvalue']
    results["SBS7b.observed"] = mutsig['SBS7b.observed']
    results["SBS7b.pvalue"] = mutsig['SBS7b.pvalue']
    results["SBS7c.observed"] = mutsig['SBS7c.observed']
    results["SBS7c.pvalue"] = mutsig['SBS7c.pvalue']
    results["SBS7d.observed"] = mutsig['SBS7d.observed']
    results["SBS7d.pvalue"] = mutsig['SBS7d.pvalue']
    results["SBS8.observed"] = mutsig['SBS8.observed']
    results["SBS8.pvalue"] = mutsig['SBS8.pvalue']
    results["SBS9.observed"] = mutsig['SBS9.observed']
    results["SBS9.pvalue"] = mutsig['SBS9.pvalue']
    results["SBS10a.observed"] = mutsig['SBS10a.observed']
    results["SBS10a.pvalue"] = mutsig['SBS10a.pvalue']
    results["SBS10b.observed"] = mutsig['SBS10b.observed']
    results["SBS10b.pvalue"] = mutsig['SBS10b.pvalue']
    results["SBS11.observed"] = mutsig['SBS11.observed']
    results["SBS11.pvalue"] = mutsig['SBS11.pvalue']
    results["SBS12.observed"] = mutsig['SBS12.observed']
    results["SBS12.pvalue"] = mutsig['SBS12.pvalue']
    results["SBS13.observed"] = mutsig['SBS13.observed']
    results["SBS13.pvalue"] = mutsig['SBS13.pvalue']
    results["SBS14.observed"] = mutsig['SBS14.observed']
    results["SBS14.pvalue"] = mutsig['SBS14.pvalue']
    results["SBS15.observed"] = mutsig['SBS15.observed']
    results["SBS15.pvalue"] = mutsig['SBS15.pvalue']
    results["SBS16.observed"] = mutsig['SBS16.observed']
    results["SBS16.pvalue"] = mutsig['SBS16.pvalue']
    results["SBS17a.observed"] = mutsig['SBS17a.observed']
    results["SBS17a.pvalue"] = mutsig['SBS17a.pvalue']
    results["SBS17b.observed"] = mutsig['SBS17b.observed']
    results["SBS17b.pvalue"] = mutsig['SBS17b.pvalue']
    results["SBS18.observed"] = mutsig['SBS18.observed']
    results["SBS18.pvalue"] = mutsig['SBS18.pvalue']
    results["SBS19.observed"] = mutsig['SBS19.observed']
    results["SBS19.pvalue"] = mutsig['SBS19.pvalue']
    results["SBS20.observed"] = mutsig['SBS20.observed']
    results["SBS20.pvalue"] = mutsig['SBS20.pvalue']
    results["SBS21.observed"] = mutsig['SBS21.observed']
    results["SBS21.pvalue"] = mutsig['SBS21.pvalue']
    results["SBS22.observed"] = mutsig['SBS22.observed']
    results["SBS22.pvalue"] = mutsig['SBS22.pvalue']
    results["SBS23.observed"] = mutsig['SBS23.observed']
    results["SBS23.pvalue"] = mutsig['SBS23.pvalue']
    results["SBS24.observed"] = mutsig['SBS24.observed']
    results["SBS24.pvalue"] = mutsig['SBS24.pvalue']
    results["SBS25.observed"] = mutsig['SBS25.observed']
    results["SBS25.pvalue"] = mutsig['SBS25.pvalue']
    results["SBS26.observed"] = mutsig['SBS26.observed']
    results["SBS26.pvalue"] = mutsig['SBS26.pvalue']
    results["SBS27.observed"] = mutsig['SBS27.observed']
    results["SBS27.pvalue"] = mutsig['SBS27.pvalue']
    results["SBS28.observed"] = mutsig['SBS28.observed']
    results["SBS28.pvalue"] = mutsig['SBS28.pvalue']
    results["SBS29.observed"] = mutsig['SBS29.observed']
    results["SBS29.pvalue"] = mutsig['SBS29.pvalue']
    results["SBS30.observed"] = mutsig['SBS30.observed']
    results["SBS30.pvalue"] = mutsig['SBS30.pvalue']
    results["SBS31.observed"] = mutsig['SBS31.observed']
    results["SBS31.pvalue"] = mutsig['SBS31.pvalue']
    results["SBS32.observed"] = mutsig['SBS32.observed']
    results["SBS32.pvalue"] = mutsig['SBS32.pvalue']
    results["SBS33.observed"] = mutsig['SBS33.observed']
    results["SBS33.pvalue"] = mutsig['SBS33.pvalue']
    results["SBS34.observed"] = mutsig['SBS34.observed']
    results["SBS34.pvalue"] = mutsig['SBS34.pvalue']
    results["SBS35.observed"] = mutsig['SBS35.observed']
    results["SBS35.pvalue"] = mutsig['SBS35.pvalue']
    results["SBS36.observed"] = mutsig['SBS36.observed']
    results["SBS36.pvalue"] = mutsig['SBS36.pvalue']
    results["SBS37.observed"] = mutsig['SBS37.observed']
    results["SBS37.pvalue"] = mutsig['SBS37.pvalue']
    results["SBS38.observed"] = mutsig['SBS38.observed']
    results["SBS38.pvalue"] = mutsig['SBS38.pvalue']
    results["SBS39.observed"] = mutsig['SBS39.observed']
    results["SBS39.pvalue"] = mutsig['SBS39.pvalue']
    results["SBS40.observed"] = mutsig['SBS40.observed']
    results["SBS40.pvalue"] = mutsig['SBS40.pvalue']
    results["SBS41.observed"] = mutsig['SBS41.observed']
    results["SBS41.pvalue"] = mutsig['SBS41.pvalue']
    results["SBS42.observed"] = mutsig['SBS42.observed']
    results["SBS42.pvalue"] = mutsig['SBS42.pvalue']
    results["SBS43.observed"] = mutsig['SBS43.observed']
    results["SBS43.pvalue"] = mutsig['SBS43.pvalue']
    results["SBS44.observed"] = mutsig['SBS44.observed']
    results["SBS44.pvalue"] = mutsig['SBS44.pvalue']
    results["SBS45.observed"] = mutsig['SBS45.observed']
    results["SBS45.pvalue"] = mutsig['SBS45.pvalue']
    results["SBS46.observed"] = mutsig['SBS46.observed']
    results["SBS46.pvalue"] = mutsig['SBS46.pvalue']
    results["SBS47.observed"] = mutsig['SBS47.observed']
    results["SBS47.pvalue"] = mutsig['SBS47.pvalue']
    results["SBS48.observed"] = mutsig['SBS48.observed']
    results["SBS48.pvalue"] = mutsig['SBS48.pvalue']
    results["SBS49.observed"] = mutsig['SBS49.observed']
    results["SBS49.pvalue"] = mutsig['SBS49.pvalue']
    results["SBS50.observed"] = mutsig['SBS50.observed']
    results["SBS50.pvalue"] = mutsig['SBS50.pvalue']
    results["SBS51.observed"] = mutsig['SBS51.observed']
    results["SBS51.pvalue"] = mutsig['SBS51.pvalue']
    results["SBS52.observed"] = mutsig['SBS52.observed']
    results["SBS52.pvalue"] = mutsig['SBS52.pvalue']
    results["SBS53.observed"] = mutsig['SBS53.observed']
    results["SBS53.pvalue"] = mutsig['SBS53.pvalue']
    results["SBS54.observed"] = mutsig['SBS54.observed']
    results["SBS54.pvalue"] = mutsig['SBS54.pvalue']
    results["SBS55.observed"] = mutsig['SBS55.observed']
    results["SBS55.pvalue"] = mutsig['SBS55.pvalue']
    results["SBS56.observed"] = mutsig['SBS56.observed']
    results["SBS56.pvalue"] = mutsig['SBS56.pvalue']
    results["SBS57.observed"] = mutsig['SBS57.observed']
    results["SBS57.pvalue"] = mutsig['SBS57.pvalue']
    results["SBS58.observed"] = mutsig['SBS58.observed']
    results["SBS58.pvalue"] = mutsig['SBS58.pvalue']
    results["SBS59.observed"] = mutsig['SBS59.observed']
    results["SBS59.pvalue"] = mutsig['SBS59.pvalue']
    results["SBS60.observed"] = mutsig['SBS60.observed']
    results["SBS60.pvalue"] = mutsig['SBS60.pvalue']
    results["SBS84.observed"] = mutsig['SBS84.observed']
    results["SBS84.pvalue"] = mutsig['SBS84.pvalue']
    results["SBS85.observed"] = mutsig['SBS85.observed']
    results["SBS85.pvalue"] = mutsig['SBS85.pvalue']

if HLAoutput is not None:
    ## parse polysovler winners.hla.txt
    winners = pd.read_csv(HLAoutput, sep="\t", header=None)
    ## create HLA columns
    results["HLA-A1"] = [winners.loc[0][1]]
    results["HLA-A2"] = [winners.loc[0][2]]
    results["HLA-B1"] = [winners.loc[1][1]]
    results["HLA-B2"] = [winners.loc[1][2]]
    results["HLA-C1"] = [winners.loc[2][1]]
    results["HLA-C2"] = [winners.loc[2][2]]

## add LOHHLA results here

## TMB calculation

if MAF_input is not None and coding_baits_BED is not None:
    ## read in MAF
    maf = pd.read_csv(MAF_input, sep="\t")
    ## prevent odd warning based on numpy type differences
    maf.Mutation_Status = maf.Mutation_Status.astype(str)
    ## subset based on only non-synonymous coding mutations
    maf = maf[maf.Mutation_Status != 'GERMLINE']
    variant_classes = ["Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation", "Frame_Shift_Ins", "Frame_Shift_Del","In_Frame_Del","In_Frame_Ins","Translation_Start_Site", "Splice_Site", "Splice_Region"]
    non_syn_mut = maf[maf.Variant_Classification.isin(variant_classes)].reset_index(drop=True)
    ## only take Chromosome, Start_Position, End_Position
    non_syn_mut_bed = non_syn_mut[["Chromosome", "Start_Position", "End_Position", "Variant_Classification", "Hugo_Symbol"]]  ## keeping the last two columns for debugging
    ## read in coding_baits_BE
    coding_baits_regions = pd.read_csv(coding_baits_BED, sep="\t", header=None) 
    coding_baits_regions.rename(columns={0:'Chromosome', 1:'Start_Position', 2:'End_Position'}, inplace=True)
    ## there is a way to do this with pandas, but others can investigate; this wrapper around bedtools is quite quic
    nonSynonymousMuts = pybedtools.BedTool.from_dataframe(non_syn_mut_bed)
    CodingRegionsBaits = pybedtools.BedTool.from_dataframe(coding_baits_regions)
    resulting_intersection = nonSynonymousMuts.intersect(CodingRegionsBaits)
    ## read into a pandas DataFrame, and count the number of rows; this is the number of mutations in non-synonymous coding mutations in canonical exons only
    mutationNum = 0
    if os.path.getsize(resulting_intersection.fn) > 0:
        resultdf = pd.read_csv(resulting_intersection.fn, sep="\t", header=None)
	mutationNum = len(resultdf.index)
    ## conditional based on BED file used:
    ## AgilentExon_51MB_b37_v3_baits.bed, total_cds_size = 30.89918
    ## IDT_Exome_v1_FP_b37_baits.bed, total_cds_size = 36.00458
    ## WGS, total_cds_size = 45.57229
    if "AgilentExon_51MB_b37" in coding_baits_BED:
    	tmb = mutationNum/30.89918
    elif "IDT_Exome_v1_FP_b37" in coding_baits_BED:
    	tmb = mutationNum/36.00458
    elif "b37_wgs_calling_regions" in coding_baits_BED:
    	tmb = mutationNum/45.57229
    else: 
    	tmb = None ## this shouldn't happen

    results["Number_of_Mutations"] = len(maf.index)
    results['TMB']=[tmb]

## write to *tsv
results.to_csv(str(sampleID + '_metadata.txt'), sep="\t", index=False)
