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

def run_pair_to_bed(bedpe,regions,match_type):
	"""
	use pybedtools to run bedtools pairtobed
	a length of 1 base must be artificially added to each end with a length of 0, otherwise no overlap
	"""
	[meta_header, bedpe_header_list] = get_bedpe_header(bedpe)
	a = bedtool_to_df(pybedtools.BedTool(bedpe), bedpe_header_list)
	if a.shape[0] > 0:
		a["END_A"] = a.apply(lambda row: row["END_A"] + 1 if row["START_A"] == row["END_A"] else row["END_A"],axis=1)
		a["END_B"] = a.apply(lambda row: row["END_B"] + 1 if row["START_B"] == row["END_B"] else row["END_B"],axis=1)
	a = pybedtools.BedTool.from_dataframe(a)
	b = pybedtools.BedTool(regions)
	result = a.pair_to_bed(b, **{'type': match_type})
	return result

def bedtool_to_df(bt,header_list):
	try:
		df = bt.to_dataframe(header=None, comment="#")
		df = df[df.columns[:len(header_list)]]
		df.columns = header_list
	except Exception as e:
		print(e)
		print("unable to convert bedtool to pandas dataframe, continuing with empty dataframe")
		df = pd.DataFrame(columns = header_list)
		
	df = df.astype({i:int for i in "START_A|END_A|START_B|END_B".split("|")})
	for j in ["#CHROM_A","CHROM_B"]:
		if df[j].dtype == 'float64':
			df = df.astype({j:int}).astype({j:str})
	return df

def get_bedpe_header(bedpe):
	with open(bedpe, 'r') as f:
		x = f.readlines()
	meta_header = [x[y] for y in range(len(x)) if x[y].startswith("##")]
	header_idx = [y for y in range(len(x)) if x[y].startswith("#CHROM")][0]
	bedpe_header_list = x[header_idx].rstrip().split("\t")
	return [meta_header, bedpe_header_list]
