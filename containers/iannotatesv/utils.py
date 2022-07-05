import argparse
import numpy as np
import pandas as pd
import pybedtools


overlap_type = {"both":"notboth",
				"notboth":"both",
				"either":"neither",
				"neither":"either"
				}

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

def parse_svtools_bedpe_file(bedpe, offset = True):
	"""
	Read bedpe file
	1. Separate components meta-data, header line, and records (main data)
	2. Coerce chromosome values to string, and position coordinates to integer
	3. if offset set to True, add +1 to END_* coordinates. This may be necessary for certain pybedtools operations.
	"""
	meta_header=""
	with open(bedpe, 'r') as f:
		main_data = False
		while main_data == False:
			x = f.readline()
			if x.startswith("##"):
				meta_header += x
			else:
				header = x
				header_list = header.strip().split("\t")
				main_data = True
		try:
			records_df = pd.read_csv(f, header=None, sep="\t" )
			records_df.columns = header_list
		except:
			records_df = pd.DataFrame(columns = header_list)

	try:
		records_df = records_df.astype({i:int for i in "START_A|END_A|START_B|END_B".split("|")})
		for j in ["#CHROM_A","CHROM_B"]:
			if records_df[j].dtype == 'float64':
				records_df = records_df.astype({j:int}).astype({j:str})
	except Exception as e:
		print(e)
		print("Unable to apply formatting to bedpe columns in pandas. Skipping...")

	try:
		if offset:
			if records_df.shape[0] > 0:
				records_df["END_A"] = records_df.apply(lambda row: row["END_A"] + 1 if row["START_A"] == row["END_A"] else row["END_A"],axis=1)
				records_df["END_B"] = records_df.apply(lambda row: row["END_B"] + 1 if row["START_B"] == row["END_B"] else row["END_B"],axis=1)
	except Exception as e:
		print(e)
		print("Unable to add +1 offset to END_A and END_B columns in pandas. There may be issues with pybedtools. Skipping...")
	return [meta_header, header_list, records_df]
