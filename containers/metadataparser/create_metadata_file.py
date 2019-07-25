#!/usr/bin/env python3


__author__  = "Evan Biederstedt"
__email__   = "biederse@mskcc.org; evan.biederstedt@gmail.com"
__version__ = "0.2.0"
__status__  = "Dev"


"""
Parse pipeline outputs into a single *tsv file, a tab-delimited *tsv of metadata

Usage: python3 create_metadata_file.py [-h] --sampleID SAMPLEID
                               [--facetsPurity_out FACETSPURITY_OUT]
                               [--facetsArmLevel FACETSARMLEVEL]
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
     --facetsArmLevel FACETSARMLEVEL
                            FACETS armLevel output, *armlevel.tsv; used to measure WGD
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

parser = argparse.ArgumentParser()
parser.add_argument('--sampleID', help = 'sample ID from channel, should always exist', required = True)
parser.add_argument('--facetsPurity_out', help = 'FACETS purity output, *_purity.out; used to measure ploidy and purity', required = False)
parser.add_argument('--facetsArmLevel', help = 'FACETS armLevel output, *armlevel.tsv; used to measure WGD', required = False)
parser.add_argument('--MSIsensor_output', help = 'MSIsensor output, *.msisensor.tsv', required = False)
parser.add_argument('--mutational_signatures_output', help = 'output *txt of mutational signatures', required = False)
parser.add_argument('--polysolver_output', help = 'output *txt of HLA genotypes, winners.hla.txt', required = False)
parser.add_argument('--MAF_input', help = 'annotated MAF, *maf', required = False)
parser.add_argument('--coding_baits_BED', help = 'BED of assay-specific baits intersecting coding regions, ensGene.all_CODING_exons.reference.bed', required = False)


## NOTE: only three asssys/BED files are acceptable for coding_baits_BED:
### --- Agilent: AgilentExon_51MB_b37_v3_baits.coding.sorted.merged.bed
### --- IDT: IDT_Exome_v1_FP_b37_baits.coding.sorted.merged.bed
### --- WGS: b37_wgs_calling_regions.v1.coding.sorted.merged.bed

args = parser.parse_args()

sampleID = args.sampleID
facetsPurityPloidy = args.facetsPurity_out
facetsArmLevel = args.facetsArmLevel
MSIoutput = args.MSIsensor_output
mutationalSignatures = args.mutational_signatures_output
HLAoutput = args.polysolver_output
MAF_input = args.MAF_input
coding_baits_BED = args.coding_baits_BED

## Current plan: put this in a pandas DataFrame and save as a *tsv file

results = pd.DataFrame()

## create sampleID column
results = results.assign(sample=[sampleID])

if facetsPurityPloidy is not None:
    ## parse purity and ploidy from FACETS *_purity.out
    purity_out = pd.read_csv(facetsPurityPloidy, sep="\t")
    ## create purity and ploidy columns
    results["purity"] = [purity_out.iloc[15].str.split(' = ')[0][1]]
    results["ploidy"] = [purity_out.iloc[16].str.split(' = ')[0][1]]


if facetsArmLevel is not None:
    ## parse WGD from facets-suite *.armlevel.tsv
    armLevel = pd.read_csv(facetsArmLevel, sep="\t")
    ## create WGD column
    if str(armLevel.WGD.unique()).replace("['", "").replace("']", "")=='WGD':
    	results["WGD_status"] = [True]
    elif str(armLevel.WGD.unique()).replace("['", "").replace("']", "")=='no WGD':
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
    results["Number_of_Mutations"] = mutsig['Number of Mutations']
    results["Signature_1"] =  mutsig['Signature.1']
    results["Signature_2"] =  mutsig['Signature.2']
    results["Signature_3"] =  mutsig['Signature.3']
    results["Signature_4"] =  mutsig['Signature.4']
    results["Signature_5"] =  mutsig['Signature.5']
    results["Signature_6"] =  mutsig['Signature.6']
    results["Signature_7"] =  mutsig['Signature.7']
    results["Signature_8"] =  mutsig['Signature.8']
    results["Signature_9"] =  mutsig['Signature.9']
    results["Signature_10"] =  mutsig['Signature.10']
    results["Signature_11"] =  mutsig['Signature.11']
    results["Signature_12"] =  mutsig['Signature.12']
    results["Signature_13"] =  mutsig['Signature.13']
    results["Signature_14"] =  mutsig['Signature.14']
    results["Signature_15"] =  mutsig['Signature.15']
    results["Signature_16"] =  mutsig['Signature.16']
    results["Signature_17"] =  mutsig['Signature.17']
    results["Signature_18"] =  mutsig['Signature.18']
    results["Signature_19"] =  mutsig['Signature.19']
    results["Signature_20"] =  mutsig['Signature.20']
    results["Signature_21"] =  mutsig['Signature.21']
    results["Signature_22"] =  mutsig['Signature.22']
    results["Signature_23"] =  mutsig['Signature.23']
    results["Signature_24"] =  mutsig['Signature.24']
    results["Signature_25"] =  mutsig['Signature.25']
    results["Signature_26"] =  mutsig['Signature.26']
    results["Signature_27"] =  mutsig['Signature.27']
    results["Signature_28"] =  mutsig['Signature.28']
    results["Signature_29"] =  mutsig['Signature.29']
    results["Signature_30"] =  mutsig['Signature.30']

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
    variant_classes = ["Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation", "Frame_Shift_Ins", "Frame_Shift_Del","In_Frame_Del","In_Frame_Ins","Translation_Start_Site", "Splice_Site"]
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
    resultdf = pd.read_csv(resulting_intersection.fn, sep="\t", header=None)
    ## conditional based on BED file used:
    ## AgilentExon_51MB_b37_v3_baits.bed, total_cds_size = 30.89918
    ## IDT_Exome_v1_FP_b37_baits.bed, total_cds_size = 36.00458
    ## WGS, total_cds_size = 45.57229
    if "AgilentExon_51MB_b37" in coding_baits_BED:
    	tmb = len(resultdf.index)/30.89918
    elif "IDT_Exome_v1_FP_b37" in coding_baits_BED:
    	tmb = len(resultdf.index)/36.00458
    elif "b37_wgs_calling_regions" in coding_baits_BED:
    	tmb = len(resultdf.index)/45.57229
    else: 
    	tmb = None ## this shouldn't happen

    results['TMB']=[tmb]

## write to *tsv
results.to_csv(str(sampleID + '_metadata.tsv'), sep="\t", index=False)

