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
		x = 0
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
	insert_mean, insert_std = estimateInsertSizeDistribution(args.bam)
	print(insert_mean)
	print(insert_std)
	cfg.set('BamParameters', 'insert_mean', str(int(insert_mean)))
	cfg.set('BamParameters', 'insert_std',  str(int(insert_std)))
	cfg.set('SVclasses', 'itrx_class',  "INTRX,TRA,BND")
	with open(cfg_out, 'w') as configfile:
		cfg.write(configfile)

	print(svclone_inputs)

def estimateInsertSizeDistribution(bamfile, alignments=10000):
	sam = pysam.AlignmentFile(bamfile)
	inserts = np.array([read.tlen for read in sam.head(alignments) if read.tlen > 0 and read.tlen < 500000 and read.is_paired and read.next_reference_id == read.reference_id ])
	print(len(inserts))
	insert_mean, insert_std = np.mean(inserts), np.std(inserts)
	return [ insert_mean, insert_std ]

def readConfigFile(filePath):
	Config = configparser.ConfigParser(delimiters=":")
	cfg_file = Config.read(filePath)
	return cfg_file

def writeConfigFile(filePath,cfg_obj):
	with open(filePath, 'w') as configfile:
		config.write(configfile,space_around_delimiters=True)
		

if __name__ == "__main__":
	main()

