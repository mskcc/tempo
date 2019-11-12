#!/usr/bin/env python
# -*- coding: utf-8 -*-
 
"""
Retrieve original variant data for merged and annotated Delly and Manta VCF.
Usage: get-sv-info.py merge.vcf delly.vcf manta.vcf output.vcf
"""

__author__  = "Philip Jonsson"
__email__   = "jonssonp@mskcc.org"
__version__ = "0.1.0"
__status__  = "Dev"

import sys, os, re
from pysam import VariantFile    # version >= 0.15.2

# Input files
merge_vcf = VariantFile(sys.argv[1], "r")
delly_vcf = VariantFile(sys.argv[2], "r")
manta_vcf = VariantFile(sys.argv[3], "r")
normal = merge_vcf.header.samples[0]
tumor = merge_vcf.header.samples[1]

# Set output file and its header
merge_vcf.header.formats.add('PR', 2, 'Integer', 'Paired-read support for reference and alternate allele, from first caller listed in ALGORITHMS')
merge_vcf.header.formats.add('SR', 2, 'Integer', 'Split-read support for reference and alternate allele, from  first caller in ALGORITHMS')
merge_vcf.header.info.add('TUMOR_PR_OTHER', 1, 'String', 'Tumor paired-read support for reference and alternate allele, from second caller last in ALGORITHMS')
merge_vcf.header.info.add('TUMOR_SR_OTHER', 1, 'String', 'Tumor split-read support for reference and alternate allele, from second caller listed in ALGORITHMS')
merge_vcf.header.info.add('NORMAL_PR_OTHER', 1, 'String', 'Normal paired-read support for reference and alternate allele, from second caller last in ALGORITHMS')
merge_vcf.header.info.add('NORMAL_SR_OTHER', 1, 'String', 'Normal split-read support for reference and alternate allele, from second caller listed in ALGORITHMS')
merge_vcf.header.info.add('BND_DEPTH', 1, 'Integer', 'Read depth at local translocation breakend')
merge_vcf.header.info.add('MATE_BND_DEPTH', 1, 'Integer', 'Read depth at remote translocation mate breakend')
merge_vcf.header.info.add('CT', 1, 'String', 'Orientation of paired-end SV')
merge_vcf.header.info.add('INSLEN', 1, 'Integer', 'Size of insertion (Delly)')
merge_vcf.header.info.add('SVINSLEN', 1, 'Integer', 'Size of insertion (Manta)')
merge_vcf.header.info.add('SVINSSEQ', 1, 'String', 'Sequence of insertion')
merge_vcf.header.info.add('HOMLEN', 1, 'Integer', 'Length of microhomology at breakpoints')
merge_vcf.header.info.add('HOMSEQ', 1, 'String', 'Sequence of microhomology at breakpoints')
vcf_out = VariantFile(sys.argv[4], "w", header = merge_vcf.header)

# Get variants
merge_vars = [var for var in merge_vcf.fetch()]
delly_vars = {var.id: var for var in delly_vcf.fetch()}
manta_vars = {var.id: var for var in manta_vcf.fetch()}

def main():

    for var in merge_vars:
        process_variant(var)
    vcf_out.close()

def process_variant(var):

    algos = var.info['ALGORITHMS']
    ids = var.info['MEMBERS']
    if all(al in algos for al in ['delly', 'manta']):
        delly_id = [x for x in ids if 'Manta' not in x][0]
        manta_id = ids[0].replace('_', ':') if 'Manta' in ids[0] else ids[1].replace('_', ':') # these are converted by svtk
    elif 'delly' in algos:
        delly_id = ids[0]
    elif 'manta' in algos:
        manta_id = ids[0].replace('_', ':')
    
    # Extract information from original Delly VCF
    # Below, set a common format for reads support from split (SR) and paired reads (PR), with reference and variant read support as tuple
    delly_info = None
    manta_info = None
    if 'delly' in algos:

        delly_var = delly_vars[delly_id]

        # Get read support
        delly_format = delly_var.format.keys()
        keep_format = ['DR', 'DV', 'RR', 'RV'] # ref and alt support for split reads (DR, DV) and paired reads (RR, RV)
        delly_tumor_gt = [(key, delly_var.samples[tumor][key]) for key in keep_format]
        delly_normal_gt = [(key, delly_var.samples[normal][key]) for key in keep_format]
        delly_tumor_new_gt = {'PR': ([x[1] for x in delly_tumor_gt if x[0] == 'RR'][0],
                                     [x[1] for x in delly_tumor_gt if x[0] == 'RV'][0]), 
                              'SR': ([x[1] for x in delly_tumor_gt if x[0] == 'DR'][0],
                                     [x[1] for x in delly_tumor_gt if x[0] == 'DV'][0])}
        delly_normal_new_gt = {'PR': ([x[1] for x in delly_normal_gt if x[0] == 'RR'][0],
                                      [x[1] for x in delly_normal_gt if x[0] == 'RV'][0]), 
                               'SR': ([x[1] for x in delly_normal_gt if x[0] == 'DR'][0],
                                      [x[1] for x in delly_normal_gt if x[0] == 'DV'][0])}                             

        # Get INFO values
        keep_info = ['INSLEN', # length of insertion
                     'HOMLEN', # length of microhomology
                     'CT'] # orientation of event
        delly_info = {key: delly_var.info[key] for key in delly_var.info.keys() if key in keep_info}
        
    # Extract information from original Manta VCF    
    if 'manta' in algos:

        manta_var = manta_vars[manta_id]
        
        # Get read support
        manta_format = manta_var.format.keys() # Manta only records ref, alt read support (SR: split reads; PR: paired reads)
        manta_tumor_gt = [(key, manta_var.samples[tumor][key]) for key in manta_format]
        manta_normal_gt = [(key, manta_var.samples[normal][key]) for key in manta_format]
        manta_tumor_new_gt = {} # unlike Delly, these values are not always present in Manta output
        manta_normal_new_gt = {}
        if 'PR' in manta_format:
            manta_tumor_new_gt.update({'PR': [x[1] for x in manta_tumor_gt if x[0] == 'PR'][0]})
            manta_normal_new_gt.update({'PR': [x[1] for x in manta_normal_gt if x[0] == 'PR'][0]})
        else:
            manta_tumor_new_gt.update({'PR': (0, 0)})
            manta_normal_new_gt.update({'PR': (0, 0)})
        if 'SR' in manta_format:
            manta_tumor_new_gt.update({'SR': [x[1] for x in manta_tumor_gt if x[0] == 'SR'][0]})
            manta_normal_new_gt.update({'SR': [x[1] for x in manta_normal_gt if x[0] == 'SR'][0]})
        else:
            manta_tumor_new_gt.update({'SR': (0, 0)})
            manta_normal_new_gt.update({'SR': (0, 0)})
        
        # Get INFO values
        keep_info = ['BND_DEPTH', # for translocations, read depth at breakpoint
                     'MATE_BND_DEPTH', # for translocations, read depth at remote breakpoint
                     'SVLEN', # length of event, negative for deletions
                     'SVINSLEN', # length of insertion
                     'SVINSSEQ', # sequence of insertion
                     'HOMLEN', # length of microhomology at breakpoints
                     'HOMSEQ'] # microhomology sequences at breakpoints
        manta_info = {key: manta_var.info[key] for key in manta_var.info.keys() if key in keep_info}

    # Set new genotype information for tumor and normal
    # var.samples[tumor].__delitem__('GT')
    # var.samples[normal].__delitem__('GT')
    tumor_split_reads = delly_tumor_new_gt['SR'] if 'delly' in algos else manta_tumor_new_gt['SR']
    normal_split_reads = delly_normal_new_gt['SR'] if 'delly' in algos else manta_normal_new_gt['SR']
    tumor_paired_reads = delly_tumor_new_gt['PR'] if 'delly' in algos else manta_tumor_new_gt['PR']
    normal_paired_reads = delly_normal_new_gt['PR'] if 'delly' in algos else manta_normal_new_gt['PR']
    
    if all(al in algos for al in ['delly', 'manta']):
        add_info = {'TUMOR_PR_OTHER': ','.join(map(str, manta_tumor_new_gt['PR'])),
                    'TUMOR_SR_OTHER': ','.join(map(str, manta_tumor_new_gt['SR'])),
                    'NORMAL_PR_OTHER': ','.join(map(str, manta_normal_new_gt['PR'])),
                    'NORMAL_SR_OTHER': ','.join(map(str, manta_normal_new_gt['SR']))}
    else:
        add_info = dict()

    var.format.clear()
    for x in tumor, normal:
        var.samples[x].keys().append('PR')
        var.samples[x].keys().append('SR')

    var.samples[tumor].__setitem__('PR', tumor_paired_reads)
    var.samples[tumor].__setitem__('SR', tumor_split_reads)
    var.samples[normal].__setitem__('PR', normal_paired_reads)
    var.samples[normal].__setitem__('SR', normal_split_reads)

    if delly_info:
        add_info.update(delly_info)
    if manta_info:
        add_info.update(manta_info)
    
    if add_info:
        for k in add_info:
            var.info.__setitem__(k, add_info[k])

    vcf_out.write(var)

if __name__ == "__main__":
    main()
