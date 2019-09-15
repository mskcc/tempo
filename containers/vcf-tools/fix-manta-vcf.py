#!/usr/bin/env python
# -*- coding: utf-8 -*-
 
"""
Fix formatting of VCF of somatic SV calls from Manta.
Usage: fix-manta-vcf.py input_filename.vcf
Output: 'input_filename.fixed.vcf'
"""

__author__  = "Philip Jonsson"
__email__   = "jonssonp@mskcc.org"
__version__ = "0.1.0"
__status__  = "Dev"

import sys, os
from pysam import VariantFile    # version >= 0.15.2

vcf_in = VariantFile(sys.argv[1], "r")
normal = vcf_in.header.samples[0]
tumor = vcf_in.header.samples[1]

# Add header tag
vcf_in.header.formats.add('GT', 1, 'String', 'Genotype')

outfile = os.path.splitext(sys.argv[1])[0] + '.fixed.vcf'
vcf_out = VariantFile(outfile, "w", header = vcf_in.header)

for var in vcf_in.fetch():

    var.samples[tumor].__setitem__('GT', (0, 1))
    var.samples[normal].__setitem__('GT', (0, 0))

    vcf_out.write(var)

vcf_out.close()
