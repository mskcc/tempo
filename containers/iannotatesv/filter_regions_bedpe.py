import argparse, sys, os
import numpy as np
import pandas as pd
import pybedtools
from utils import *

def usage():
	parser = argparse.ArgumentParser()
	parser.add_argument('--blacklist-regions', dest="regions", help = 'regions bed/bedpe file', required = True)
	parser.add_argument('--bedpe', help = 'bedpe file', required = True)
	parser.add_argument('--tag', help = 'string to update FILTER column', required = True)
	parser.add_argument('--output',dest="outfile",help = 'output file' , required = True)
	parser.add_argument('--match-type',dest="type",help = 'both/either/notboth/neither' , required = True)
	parser.add_argument('--ignore-strand',dest="ignore_strand", default=False, action="store_true", help = 'Use flag to ignore strand when running bedtools. Default False' )
	return parser.parse_args()

def determine_regions_filetype(filepath):
	"""
	Determine if the regions file is bed or bedpe, and raise error if neither
	"""
	if not any([filepath.endswith(i) for i in "bed|bedpe|bed.gz|bedpe.gz".split("|")]):
		raise ValueError("--bedpe requires bed, bedpe, bed.gz or bedpe.gz file.")
	elif any([filepath.endswith(i) for i in "bed|bed.gz".split("|")]):
		return True
	else: return False

def validate_bedpe_input(filepath):
	"""
	Raise error if filetype is not bedpe
	"""
	if not any([filepath.endswith(i) for i in "bedpe|bedpe.gz".split("|")]):
		raise ValueError("--bedpe requires bedpe or bedpe.gz file.")

def validate_match_type(type):
	"""
	Raise error if type is not acceptable for use in pybedtools pair_to_pair or pair_to_bed
	"""
	if type not in overlap_type.keys():
		raise ValueError( "--match-type must be {}.".format(" or ".join(overlap_type.keys())) )

def main():
	"""
	1. validate inputs
	2. read input beds/bedpes as pybedtools.BedTool objects
	3. extract variant IDs from overlapping regions
	4. update the FILTER tag for identified variants
	"""
	args = usage()

	# parse inputs
	try:
		validate_bedpe_input(args.bedpe)
		[meta_header, bedpe_header_list, bedpe_df] = parse_svtools_bedpe_file(args.bedpe)
	except Exception as e:
		print(e)
		sys.exit("Unable to parse --bedpe")

	try:
		is_bed = determine_regions_filetype(args.regions)
		regions_bt = pybedtools.BedTool(args.regions)
	except Exception as e:
		print(e)
		sys.exit("Unable to parse --blacklist-regions")

	# annotate bedpe
	try:
		bedpe_bt = pybedtools.BedTool.from_dataframe(bedpe_df)
		ids_df = find_overlapped_ids(bedpe_bt, regions_bt, args.type, args.ignore_strand, is_bed)
		bedpe_df = add_filter_by_id(bedpe_df,ids_df,args.tag)
	except Exception as e:
		print(e)
		sys.exit()

	# write result
	with open(args.outfile, "w") as fw:
		fw.write("".join(meta_header))
	bedpe_df.to_csv(args.outfile, header=True, index=False, sep="\t", mode="a")

def find_overlapped_ids(bedpe_bt,regions_bt,match_type,ignore_strand,is_bed):
	# determine validity of match_type
	validate_match_type(match_type)

	# run pair_to_bed if regions file is bed, pair_to_pair if regions file is bedpe
	if is_bed:
		intersect = run_pair_to_bed(bedpe_bt,regions_bt,match_type)
	else:
		intersect = run_pair_to_pair(bedpe_bt,regions_bt,match_type, ignore_strand=ignore_strand)

	# extract ids
	try:
		ids_df = intersect.to_dataframe(header=None)[6].drop_duplicates()
	except:
		ids_df = pd.DataFrame(columns=["name"])
	return ids_df

def add_filter_by_id(bedpe_df,ids_df,tag):
	# annotate bedpe based on matching IDs
	ids_df.columns = ["ID"]
	ids_df["FILTER_NEW"] = tag
	filtered_bedpe_df = bedpe_df.merge(ids_df, on="ID", how="left")
	filtered_bedpe_df = filtered_bedpe_df.apply(lambda x: update_filter(x,x["FILTER_NEW"]) if not pd.isna(x["FILTER_NEW"]) else x, axis=1)
	filtered_bedpe_df = filtered_bedpe_df.drop(['FILTER_NEW'], axis=1)

	return filtered_bedpe_dfs

if __name__ == "__main__":
	main()
