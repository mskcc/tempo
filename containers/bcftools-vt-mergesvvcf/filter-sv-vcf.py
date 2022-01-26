#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import print_function
import argparse

"""
Perform custom filtering/annotation of variants in VCF file, after merging SVs with mergesvvcf
Usage: filter-sv-vcf.py -h 

Justification:
Mergesvvcf will retain all variants that are present in a minimum of callers, but 
does not have special handling for filters from the individual callers. For example:
1       88333444        .       N       N]1:88335187]   .       manta_MinSomaticScore   END=88335187;Callers=brass,manta,svaba,delly;NumCallers=4 ...
This is called as FAILED even though 1 caller filtered it and it is pass in the other 3. 
Similarly:
7       48677419        .       N       N[7:48677548[   .       manta_MinSomaticScore   END=48677548;Callers=manta,svaba;NumCallers=2 ...
This is retained because it is present in two callers, but it is only PASS in one caller.

We need to mark calls as PASS if they are PASS in a minimum number of individual callers,
while retaining their original filter information as annotation.
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
	parser.add_argument('--min',default = 1, help = 'minimum number of PASS callers needed to PASS final variant' , required = True)
	return parser.parse_args()

def main():
	args = usage()
	filter_by_pass_callers(args.input, args.output, args.min)

def filter_by_pass_callers(input,output,min_pass):
	vcf_in = VariantFile(input, "r")
	vcf_in.header.info.add(
		"NumCallersPass", 1, "Integer", "Number of callers that made this call without filters"
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
		if num_callers_pass > min_pass:
			vcf_rec.filter.add("PASS")

		# change TRA to BND for svtools
		if vcf_rec.info["SVTYPE"] == "TRA":
			 vcf_rec.info.__setitem__("SVTYPE","BND")
		
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