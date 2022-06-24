import argparse
import numpy as np
import pandas as pd
import pybedtools

def add_tag(val,tag ):
	"""
	update FILTER value 
	"""
	newTag = tag
	if val not in [".",None,"PASS",np.nan]:
		newTag = "{};{}".format(val,tag)
	return newTag

def update_filter(row, tag):
	"""
	return row with updated FILTER value 
	"""
	row["FILTER"] = add_tag(row["FILTER"], tag)
	return row

def run_pair_to_bed(bedpe,bed,match_type):
	"""
	use pybedtools to run bedtools pairtobed
	a length of 1 base must be artificially added to each end with a length of 0, otherwise no overlap
	"""
	result = bedpe.pair_to_bed(bed, **{'type': match_type})
	return result

def run_pair_to_pair(bedpe,filter_bedpe,match_type,ignore_strand=False):
	"""
	use pybedtools to run bedtools pairtopair
	a length of 1 base must be artificially added to each end with a length of 0, otherwise no overlap
	"""
	result = bedpe.pair_to_pair(filter_bedpe, **{'type': match_type,'is':ignore_strand})
	return result


def bedtool_to_df(bt,header_list):
	try:
		df = bt.to_dataframe(header=None)
		df = df[df.columns[:len(header_list)]]
		df.columns = header_list
	except Exception as e:
		print(e)
		print("unable to convert bedtool to pandas dataframe, continuing with empty dataframe")
		df = pd.DataFrame(columns = header_list)
	try:
		df = df.astype({i:int for i in "START_A|END_A|START_B|END_B".split("|")})
		for j in ["#CHROM_A","CHROM_B"]:
			if df[j].dtype == 'float64':
				df = df.astype({j:int}).astype({j:str})
	except Exception as e:
		print(e)
		print("Unable to apply formatting to bedpe columns using pandas, when converting pybedtool.Bedtool to pandas.DataFrame. Skipping...")
	return df.drop_duplicates()

def parse_svtools_bedpe_file(bedpe):
	with open(bedpe, 'r') as f:
		x = f.readlines()
	meta_header = [x[y] for y in range(len(x)) if x[y].startswith("##")]
	header_idx = [y for y in range(len(x)) if x[y].startswith("#CHROM")][0]
	bedpe_header_list = x[header_idx].rstrip().split("\t")
	bedpe_df = pd.read_csv(bedpe, skiprows=header_idx,header=0, sep="\t")
	try:
		bedpe_df = bedpe_df.astype({i:int for i in "START_A|END_A|START_B|END_B".split("|")})
		for j in ["#CHROM_A","CHROM_B"]:
			if bedpe_df[j].dtype == 'float64':
				bedpe_df = bedpe_df.astype({j:int}).astype({j:str})
	except Exception as e:
		print(e)
		print("Unable to apply formatting to bedpe columns in pandas. Skipping...")

	try:
		if bedpe_df.shape[0] > 0:
			bedpe_df["END_A"] = bedpe_df.apply(lambda row: row["END_A"] + 1 if row["START_A"] == row["END_A"] else row["END_A"],axis=1)
			bedpe_df["END_B"] = bedpe_df.apply(lambda row: row["END_B"] + 1 if row["START_B"] == row["END_B"] else row["END_B"],axis=1)
	except Exception as e:
		print(e)
		print("Unable to add +1 offset to END_A and END_B columns in pandas. There may be issues with pybedtools. Skipping...")
	return [meta_header, bedpe_header_list,bedpe_df]
