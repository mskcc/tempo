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

	[meta_header, bedpe_header_list] = get_bedpe_header(args.bedpe)
	filtered_bed_df = filter_bedpe(args.bedpe,bedpe_header_list,args.regions,args.tag,args.type)

	with open(args.outfile, "w") as fw:
		fw.write("".join(meta_header))
	print(filtered_bed_df.head())
	filtered_bed_df.to_csv(args.outfile, header=True, index=False, sep="\t", mode="a")

def filter_bedpe(bedpe,bedpe_header_list,regions,tag_str,match_type="either"):
	filter_col = bedpe_header_list.index("FILTER")

	intersect = bedtool_to_df(run_pair_to_bed(bedpe,regions,match_type), bedpe_header_list)
	outersect = bedtool_to_df(run_pair_to_bed(bedpe,regions, "notboth" if match_type == "both" else "neither"), bedpe_header_list)
	
	intersect["FILTER"] = intersect["FILTER"].apply(lambda x: add_tag(x,tag_str))
	
	return pd.concat([intersect,outersect])

if __name__ == "__main__":
	main()
