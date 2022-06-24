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
	parser.add_argument('--match-type',dest="type",help = 'both/either' , required = True)
	parser.add_argument('--ignore-strand',dest="ignore_strand", default=False, action="store_true", help = 'Use flag to ignore strand when running bedtools. Default False' )
	return parser.parse_args()


def main():
	args = usage()

	if args.type not in ["both","either"]:
		sys.exit("--match-type must be both or either")

	match_type = args.type	
	# determine if regions file is bed or bedpe
	filename=os.path.basename(args.regions)
	is_bed = False
	if filename.endswith("bed") or filename.endswith("bed.gz"):
		is_bed = True
	if filename.endswith("bedpe") or filename.endswith("bedpe.gz"):
		pass
	else:
		sys.exit("--blacklist-regions must have .bed/.bed.gz/.bedpe/.bedpe.gz extension")

	[meta_header, bedpe_header_list, bedpe_df] = parse_svtools_bedpe_file(args.bedpe)
	bedpe_bt = pybedtools.BedTool.from_dataframe(bedpe_df)
	regions_bt = pybedtools.BedTool(args.regions)
	
	if is_bed:
		intersect = run_pair_to_bed(bedpe_bt,regions_bt,match_type)
	else:
		intersect = run_pair_to_pair(bedpe_bt,regions_bt,match_type, ignore_strand=args.ignore_strand)

	intersect_df = bedtool_to_df(intersect,bedpe_header_list)[["ID"]]
	intersect_df["FILTER_NEW"] = args.tag
	filtered_bedpe_df = bedpe_df.merge(intersect_df, on="ID", how="left")
	filtered_bedpe_df = filtered_bedpe_df.apply(lambda x: update_filter(x,x["FILTER_NEW"]) if not pd.isna(x["FILTER_NEW"]) else x, axis=1)
	filtered_bed_df = filtered_bedpe_df.drop(['FILTER_NEW'], axis=1)

	with open(args.outfile, "w") as fw:
		fw.write("".join(meta_header))
	filtered_bed_df.to_csv(args.outfile, header=True, index=False, sep="\t", mode="a")

if __name__ == "__main__":
	main()
