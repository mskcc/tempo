#!/usr/bin/env python

"""
Prepare input files (maf, bedpe, purity/ploidy file,
and config) so that it can be used in SVclone.
 * Maf file: should not contain filtered variants. The
 columns Chromosome, vcf_pos, t_ref_count, t_alt_count are
 retained.
 * Bedpe file: should not contain filtered variants. The
 columns #CHROM_A, START_A, STRAND_A, CHROM_B, START_B,
 STRAND_B and TYPE are retained.
 * Purity/ploidy file: can be produced either by ASCAT or by
 generate_samplestatistics.R script in the DoFacets process.
 * Config file: Template config file which is modified before
 running SVclone.
Usage: prepare_svclone_inputs.py -h
"""

__author__  = "Anne Marie Noronha"
__email__   = "noronhaa@mskcc.org"
__version__ = "0.0.1"
__status__  = "Dev"

import configparser
import argparse, os
import pysam, pandas as pd, numpy as np

def usage():
	parser = argparse.ArgumentParser(description="Prepare inputs for SVClone")
	parser.add_argument("--maf", help="Filtered maf file")
	parser.add_argument("--bedpe", help="Filtered SV calls in bedpe format")
	parser.add_argument("--sampleid", metavar="s_C_XXXXXX_XNNN_d__s_C_XXXXXX_XNNN_d")
	parser.add_argument("--purity_ploidy", help="purity and ploidy file from Facets process")
	parser.add_argument("--out_dir",default="svclone_inputs", help="Folder to write all svclone inputs")
	parser.add_argument("--cfg_template",help="Template config for editing")
	parser.add_argument("--bam",help="Tumor Bam file")
	parser.add_argument("--cnv",help="CNV file")

	return parser.parse_args()

def main():
	args = usage()
	
	if os.path.isdir(args.out_dir):
		pass
	else:
		os.mkdir(args.out_dir)
	mut_out = os.path.join(args.out_dir,"callstats.txt")
	sv_out = os.path.join(args.out_dir,"simple.sv.txt")
	purity_ploidy_out = os.path.join(args.out_dir,"svclone_ploidy.txt")
	cfg_out = os.path.join(args.out_dir,"svclone_config.ini")

	svclone_inputs = {"mut":mut_out,"sv":sv_out, "purity_ploidy":purity_ploidy_out, "cfg":cfg_out}

	## read in maf and output reformatted SNPs
	mut = pd.read_csv(args.maf,sep="\t",header=0)
	mut["judgement"] = "KEEP"
	mut = mut["Chromosome,vcf_pos,t_ref_count,t_alt_count,judgement".split(",")]
	mut.columns = "contig,position,t_ref_sum,t_alt_sum,judgement".split(",")
	mut.to_csv(mut_out, sep="\t", header=True, index=False)

	## read in bedpe and output reformatted SVs
	with open(args.bedpe,"r") as f:
		in_meta=True
		line_cursor = 0
		while in_meta:
			bedpe_line = f.readline()
			if bedpe_line.startswith("#CHROM"):
				in_meta=False
				f.seek(line_cursor)
				sv = pd.read_csv(f, header=0, sep="\t")
			line_cursor = f.tell()

	sv = sv["#CHROM_A,START_A,STRAND_A,CHROM_B,START_B,STRAND_B,TYPE".split(",")]
	sv.columns = "chr1,pos1,dir1,chr2,pos2,dir2,classification".split(",")
	sv.to_csv(sv_out,header=True, sep="\t", index=False )

	## read in ploidy file and output reformatted information
	purity_ploidy = pd.read_csv(args.purity_ploidy, header=None, sep=" ", index_col=0, usecols=[0, 1]).T
	purity_ploidy["sample"] = args.sampleid
	purity_ploidy = purity_ploidy[["sample","rho","Ploidy"]]
	purity_ploidy.columns = "sample,purity,ploidy".split(",")
	purity_ploidy.to_csv(purity_ploidy_out, index=False, header=True, sep="\t", )

	## read in config and adjust insert_mean, insert_std
	## also adjust SVClass descriptor terms.
	cfg = configparser.ConfigParser()
	cfg.read(args.cfg_template)
	max_cn = max(pd.read_csv(args.cnv, sep=",",header=None).iloc[:, 6].tolist())
	insert_mean, insert_std = estimateInsertSizeDistribution(args.bam)
	read_len = estimateReadLen(args.bam)
	mean_cov = estimateCoverage(args.bam)
	cfg.set('BamParameters', 'read_len', str(int(read_len)))
	cfg.set('BamParameters', 'insert_mean', str(int(insert_mean)))
	cfg.set('BamParameters', 'insert_std',  str(int(insert_std)))
	cfg.set('BamParameters', 'max_cn',  str(int(max_cn)))
	cfg.set('BamParameters', 'mean_cov',  str(int(mean_cov)))
	cfg.set('SVclasses', 'itrx_class',  "INTRX,TRA,BND")
	cfg.set('SVannotateParameters','sv_class_field','classification')
	with open(cfg_out, 'w') as configfile:
		cfg.write(configfile)

	print(svclone_inputs)

def estimateInsertSizeDistribution(bamfile, alignments=10000):
	sam = pysam.AlignmentFile(bamfile)
	inserts = np.array([read.tlen for read in sam.head(alignments) if read.tlen > 0 and read.tlen < 500000 and read.is_paired and read.next_reference_id == read.reference_id ])
	insert_mean, insert_std = np.mean(inserts), np.std(inserts)
	return [ insert_mean, insert_std ]

def readConfigFile(filePath):
	config = configparser.ConfigParser(delimiters=":")
	config.read(filePath)
	return config

def writeConfigFile(outFilePath,cfg_obj):
	with open(outFilePath, 'w') as configfile:
		cfg_obj.write(configfile,space_around_delimiters=True)
		
def estimateCoverage(bamfile):
	"""
	estimate coverage in bam file
	assumes uniform coverage
	"""
	sam = pysam.AlignmentFile(bamfile)
	cov = list()
	interval_len = 1000000
	for i in sam.references:
		if i.startswith("GL") or i.startswith("NC"): continue
		if i in ["hs37d5","X","Y","chrX","chrY"]: continue
		print(i)
		for j in range(1,int(sam.get_reference_length(i)/interval_len)):
			try:
				x = sam.count_coverage(i,start=j*interval_len, stop=(j*interval_len) + 50)
				cov += list(np.add(np.add(x[0],x[1]),np.add(x[2],x[3])))
			except ValueError as e:
				print(e)
				pass
	print(len(cov))
	return int(sum(cov)/len(cov))

def estimateReadLen(bamfile,alignments=1000):
	sam = pysam.AlignmentFile(bamfile)
	sizes = [read.rlen for read in sam.head(alignments)]
	sizes_dict = dict()
	for i in sizes:
		if not i in sizes_dict: sizes_dict[i] = 0
		sizes_dict[i] += 1
	read_len = max(sizes_dict, key=sizes_dict.get)
	if sizes_dict[read_len] < alignments * .5:
		read_len = int(sum(sizes)/len(sizes))
	return read_len

if __name__ == "__main__":
	main()