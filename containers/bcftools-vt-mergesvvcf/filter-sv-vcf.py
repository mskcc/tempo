#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import print_function
import argparse

"""
Perform custom filtering/annotation of variants in 
VCF file, after merging SVs with mergesvvcf. Variants
with a minimum of PASSing callers should be filtered 
as PASS with all non-PASS filters recorded in INFO.
Additionally, TRA is converted to BND because svtools
will not correctly handle TRA.
Usage: filter-sv-vcf.py -h 
"""

__author__  = "Anne Marie Noronha"
__email__   = "noronhaa@mskcc.org"
__version__ = "0.0.1"
__status__  = "Dev"

import sys, os
from pysam import VariantFile    # version >= 0.15.2
from itertools import groupby


def get_flag_sources(filters, callers):
	parsed_filters = [(i,j[len(i + "_"):]) for i in callers for j in filters if j.startswith(i + "_")]
	dict_1=dict()
	for caller,val in parsed_filters:
		dict_1.setdefault(caller, []).append(val)
	return dict_1

def append_file_to_file(_from, _to):
	block_size = 1024*1024
	with open(_from, "rb") as infile, open(_to,"ab") as outfile:
		while True:
			input_block = infile.read(block_size)
			if not input_block:
				break
			outfile.write(input_block)
	infile.close()
	outfile.close()

def usage():
	parser = argparse.ArgumentParser()
	parser.add_argument('--input', help = 'input file', required = True)
	parser.add_argument('--output', help = 'output file', required = True)
	parser.add_argument('--min',type=int, default = 1, help = 'minimum number of PASS callers needed to PASS final variant' , required = True)
	return parser.parse_args()

def main():
	args = usage()
	filter_by_pass_callers(args.input, args.output, args.min)

def filter_by_pass_callers(input,output,min_pass):
	vcf_in = VariantFile(input, "r")
	vcf_in.header.info.add(
		"NumCallersPass", 1, "Integer", "Number of callers that made this call without filters"
    )
	vcf_in.header.filters.add(
		"minPassFilter", None,None, "Number of callers insufficient"
    )

	all_callers = set()

	out_vcf_recs = open("vcf.records.tmp", 'w') 
	counter = 0
	for vcf_rec in vcf_in.fetch():
		counter +=1
		## Variant info
		info = vcf_rec.info.keys()
		filter = vcf_rec.filter.keys()
		new_flags = []
		callers = list(vcf_rec.info["Callers"])
		num_callers = vcf_rec.info['NumCallers']
		for i in set(callers) - all_callers:
			vcf_in.header.info.add(
				"{}_filters".format(i),
				1,
				'String', 
				'Filter values from {} caller'.format(i)
			)
		all_callers.update(set(callers))

		#parse filters
		flag_sources = get_flag_sources(filter,callers)
		num_callers_pass = num_callers - len(flag_sources)
		#annotate with number of passing callers and filter values
		vcf_rec.info.__setitem__('NumCallersPass',num_callers_pass)
		for i in flag_sources:
			vcf_rec.info.__setitem__(i + "_filters",",".join(flag_sources[i]))
			
		#adjust filter based on passing callers
		if num_callers_pass >= min_pass:
			vcf_rec.filter.add("PASS")
		else:
			vcf_rec.filter.add("minPassFilter")

		# change TRA to BND for svtools
		if vcf_rec.info["SVTYPE"] == "TRA":
			vcf_rec.info.__setitem__("SVTYPE","BND")

		# The following lines are intended to stop the loss of END= in the output. 
		# Later versions of pysam will have fixed this, this code should fixed when pysam upgraded.
		if not vcf_rec.stop or vcf_rec.stop < 1:
			#vcf_rec.stop = vcf_rec.pos + 1
			continue
		if vcf_rec.chrom == vcf_rec.info["CHR2"] and abs(vcf_rec.start - vcf_rec.stop) <= 1 :
			#vcf_rec.stop = vcf_rec.pos + 1
			continue

		
		# print result to tmp file
		print(vcf_rec, end="", file=out_vcf_recs)
	out_vcf_recs.close()

	#write header
	out_vcf_header = open("vcf.header.tmp", 'w') 
	print(vcf_in.header, end="", file=out_vcf_header)
	out_vcf_header.close()
	
	#paste together header and records in order
	if os.path.exists(output):
		os.remove(output)
	for infile in ["vcf.header.tmp","vcf.records.tmp"]:
		append_file_to_file(infile, output)

if __name__ == "__main__":
	main()
