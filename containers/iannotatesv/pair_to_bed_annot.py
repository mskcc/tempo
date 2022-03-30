import argparse, sys
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
	parser.add_argument('--match-type',dest="type",help = 'both/either' , required = True)

	return parser.parse_args()

def main():
	args = usage()

	if args.type not in ["both","either"]:
		sys.exit("--match-type must be both or either")

	[meta_header, bedpe_header_list, bedpe_df] = parse_svtools_bedpe_file(args.bedpe)
	bedpe_bt = pybedtools.BedTool.from_dataframe(bedpe_df)
	regions_bt = pybedtools.BedTool(args.regions)
	intersect,outersect = filter_bedpe(bedpe_bt,regions_bt,args.type)
	intersect_df = bedtool_to_df(intersect,bedpe_header_list)
	intersect_df["FILTER"] = intersect_df["FILTER"].apply(lambda x: add_tag(x,args.tag))
	outersect_df = bedtool_to_df(outersect,bedpe_header_list)
	filtered_bed_df = pd.concat([intersect_df,outersect_df])

	with open(args.outfile, "w") as fw:
		fw.write("".join(meta_header))
	filtered_bed_df.to_csv(args.outfile, header=True, index=False, sep="\t", mode="a")

def filter_bedpe(bedpe_bt,regions_bt,match_type="either"):
	outersect_match_type = "notboth" if match_type == "both" else "neither"
	intersect = run_pair_to_bed(bedpe_bt,regions_bt,match_type)
	outersect = run_pair_to_bed(bedpe_bt,regions_bt, outersect_match_type)
	#intersect["FILTER"] = intersect["FILTER"].apply(lambda x: add_tag(x,tag_str))
	
	return [intersect,outersect]
	#return pd.concat([intersect,outersect])

if __name__ == "__main__":
	main()
