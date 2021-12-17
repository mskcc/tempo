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
		newTag = "{},{}".format(val,tag)
	return newTag

def update_filter(row, tag):
	"""
	return row with updated FILTER value 
	"""
	row["FILTER"] = add_tag(row["FILTER"], tag)
	return row


def run_pair_to_bed(bedpe,regions,match_type):
	a = pybedtools.BedTool(bedpe)
	b = pybedtools.BedTool(regions)
	result = a.pair_to_bed(b, **{'type': match_type})
	return result
	#in_type = match_type
	#out_type = "notboth" if match_type == "both" else "neither"
	#intersect = a.pair_to_bed(b, **{'type': in_type})
	#outersect = a.pair_to_bed(b,**{'type': out_type})
	#return [intersect, outersect]

def bedtool_to_df(bt,header_list):
	df = bt.to_dataframe(header=None, comment="#")
	if df.shape[0] > 0:
		df = df[df.columns[:len(header_list)]]
		df.columns = header_list
	else:
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
