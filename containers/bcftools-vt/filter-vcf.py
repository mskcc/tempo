#!/usr/bin/env python
# -*- coding: utf-8 -*-
 
"""
Perform custom annotation of variants in VCF file, based on MuTect2 and Strelka2 variant calls and custom pre-processing.
"""

__author__  = "Philip Jonsson"
__email__   = "jonssonp@mskcc.org"
__version__ = "0.1.0"
__status__  = "Dev"

import sys, os
from pysam import VariantFile
from itertools import groupby

vcf_in = VariantFile(sys.argv[1], "r")
normal = vcf_in.header.samples[0]
tumor = vcf_in.header.samples[1]

# Add new headers
# vcf_in.header.filters.add('multiallelic', None, None, 'Multiple alleles at same locus') # Note that MuTect2 already has a FILTER tag for this, we're hijacking that
vcf_in.header.filters.add('multiallelic2', None, None, 'Multiple alleles at same locus') # Note that MuTect2 already has a FILTER tag for this, we're hijacking that
vcf_in.header.filters.add('part_of_mnv', None, None, 'Variant is part of previous variant')
# vcf_in.header.filters.add('short_repeat', None, None, 'Variant part of a short repeat')
vcf_in.header.filters.add('strand_bias', None, None, 'Variant part of a short repeat') # Note that MuTect2 already has a FILTER tag for this, we're hijacking that
# vcf_in.header.filters.add('caller_conflict', None, None, 'MuTect2 and Strelka2 provides conflicting FILTER flags for this variant')
# vcf_in.header.info.add('Caller', 1, 'String', 'Comman-separated list of variant callers that called variant')
vcf_in.header.info.add('Ref_Tri', 1, 'String', 'Normalize trinucleotide context of SNVs')


outfile = os.path.splitext(sys.argv[1])[0] + '.filter.vcf'
vcf_out = VariantFile(outfile, "w", header = vcf_in.header)
prev_var = None
for var in vcf_in.fetch():

    # Variant info
    info = var.info.keys()
    pos = var.pos
    print(var.pos)
    ref = var.ref
    alt = var.alts[0]
    filter = var.filter.keys()
    new_flags = []

    # Add caller tag
    # caller = []
    # if 'MuTect2' in info:
    #     caller.append('MuTect2')
    # if 'Strelka2' in info:
    #     caller.append('Strelka2')
    # var.info.__setitem__('Caller', ','.join(caller))

    # Add normalized reference trinucleotide context
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

    if len(ref) == len(alt) == 1:
        ref_tri = var.info['FLANKSEQ'][9:14].replace('[', '').replace(']', '')
        if ref not in ['C', 'T']:
            ref_tri = ''.join([complement[nt] for nt in ref_tri])[::-1]
        var.info.__setitem__('Ref_Tri', ref_tri)

    # Check for variant substitutions that are part of the previous variant
    # These are most likely (but not always?) Strelka2 calls that MuTect2 called as MNVs
    if prev_var is not None:
        prev_end = prev_var.pos + len(prev_var.alts[0]) - 1
        if pos == prev_end:
            ref_match = ref[0] == prev_var.ref[-1]
            alt_match = alt[0] == prev_var.alts[0][-1]
        
            if ref_match and alt_match:
                new_flags.append('part_of_mnv')
        
    # Check for likely artifacts in short repeats
    # left_flank = var.info['FLANKSEQ'].split('[')[0][::-1]
    # right_flank = var.info['FLANKSEQ'].split(']')[1]
    # if len(alt) > len(ref): # insertion
    #     alt_repeat = alt[1:]
    # elif len(alt) < len(ref): # deletion
    #     alt_repeat = ref[1:]
    # else:
    #     alt_repeat = alt

    # right_flank = [right_flank[i:i+len(alt_repeat)] for i in range(0, len(right_flank), len(alt_repeat))]
    # left_flank = [left_flank[i:i+len(alt_repeat)] for i in range(0, len(left_flank), len(alt_repeat))]
    # right_flank = [(key, len(list(group))) for key, group in groupby(right_flank)]
    # left_flank = [(key, len(list(group))) for key, group in groupby(left_flank)]
    # rep_length = len(alt_repeat) * sum([bps[1] for bps in [left_flank[0], right_flank[0]] if bps[0] == alt_repeat])
    
    # if rep_length > 5:
    #     new_flags.append('short_repeat') 

    # Check for multiallelic Strelka2 calls
    # Strelka2 does not produce multiallelic calls by itself but provides info on all alternate bases observed
    if "Strelka2" in info and len(ref) == len(alt) == 1:
        alleles = { # what about indels from Strelka2?
            'A': var.samples[tumor]['AU'],
            'C': var.samples[tumor]['CU'],
            'G': var.samples[tumor]['GU'],
            'T': var.samples[tumor]['TU']
        }        
        tier1_other = sum([alleles[a][0] for a in alleles if a not in [alt, ref]])
        tier2_other = sum([alleles[a][1] for a in alleles if a not in [alt, ref]])
        alt_reads = alleles[alt][0]

        if (tier1_other + tier2_other) >= alt_reads * 0.5:
            new_flags.append('multiallelic2')
    
    # Add an additional strand-bias filer
    # Only MuTect2 provides sufficient variant information for this
    if "MuTect2" in info:
        t_fw = var.samples[tumor]['F1R2']
        t_rev = var.samples[tumor]['F2R1']
        n_fw = var.samples[tumor]['F1R2']
        n_rev = var.samples[tumor]['F2R1']

        if t_fw[1] == 0 or t_rev[1] == 0: # if all support reads come from one read-pair orientation
            if t_fw[0] > 10 and t_rev[0] > 10 or n_fw[0] > 10 and n_rev[0] > 10:
                new_flags.append('strand_bias')

    # Conflicting MuTect2 and Strelka2 filters
    # if "PASS" in filter and "Strelka2FAIL" in info:
    #     new_flags.append('caller_conflict')

    if len(new_flags) > 0:
        if "PASS" in var.filter.keys():
            var.filter.clear()
        for flag in new_flags:
            var.filter.add(flag)    

    prev_var = var

    # Write to output
    vcf_out.write(var)

vcf_out.close()
