import argparse, sys, os
import numpy as np
import pandas as pd
import pybedtools
from utils import *

def usage():
	parser = argparse.ArgumentParser()
	parser.add_argument('--exon-junct', dest="regions", help = 'regions with exon junctions', required = True)
	parser.add_argument('--bedpe', help = 'bedpe file', required = True)
	parser.add_argument('--out',dest="outfile",help = 'output file' , required = True)
	parser.add_argument('--out-bedpe',dest="outbedpe",help = 'output bedpe' , required = True)
	parser.add_argument('--intermediate', action="store_true", help = 'print intermediate file for advanced analysis')
	return parser.parse_args()

def get_artefact_calls(df):
	"""
	Find genes with at least 2 exon-exon junctions
	remove any variant that doesn't involve two exons
	remove any variant that doesn't have "donor" and "acceptor" in the correct orientation
	Only genes with at least two splice variants will be considered. 
	"""
	try: 
		df = df.groupby(["ID","gene"]).filter(lambda x: len(set(x['exon_number'])) == 2 and len(set(x['splice_orientation'])) == 2)
		df = df.groupby(["ID","gene"]).filter(lambda x: sorted(zip(x["exon_number"],x["splice_orientation"]), key=lambda t: t[0] )[0][1] == "donor")
		df = df.groupby(["gene"]).filter(lambda x: len(set(x['ID'])) > 1)
		df = df[["ID","gene"]].drop_duplicates().rename(columns={"gene":"POTENTIAL_CDNA_CONTAMINATION"})
		return df.groupby(['ID'], as_index = False).agg({'POTENTIAL_CDNA_CONTAMINATION': ';'.join})
	except:
		return pd.DataFrame(columns=["ID","POTENTIAL_CDNA_CONTAMINATION"])
	#df = df.drop_duplicates().groupby(["var_id","gene"], as_index = False).agg({'exon': lambda x: len(set(x)) })
	#artefact_genes = df[df.exon > 3].gene.tolist()
	#return df[df.gene.isin(artefact_genes)]["var_id"].tolist()

def prep_intersection(intersect=pd.DataFrame(),bedpe_header_list=[]):
	attr_list = sorted("exon_number|splice_orientation|gene|ccds_id".split("|"))
	try:
		intersect.columns = bedpe_header_list + "#chrom|start|end|attr|score|cds_strand".split("|")
		attr_list = sorted("exon_number|splice_orientation|gene|ccds_id".split("|"))
		intersect["attr"] = intersect.apply(lambda x: {j.split(":")[0]:j.split(":")[1] for j in x["attr"].split("|")}, axis=1)
		intersect["attr"] = intersect.apply(lambda x: [ x["attr"].get(i,None) for i in attr_list  ],axis=1)
		intersect = intersect[~intersect.apply(lambda x: None in x["attr"], axis=1)]
		intersect[attr_list] = intersect.attr.apply(lambda x: pd.Series(x))
		intersect = intersect[intersect["TYPE"]=="DEL"]
		intersect = intersect[["ID"] + attr_list].drop_duplicates(keep="first")
		return intersect
	except:
		return pd.DataFrame(columns=["ID"] +attr_list)



def main():
	args = usage()

	[meta_header, bedpe_header_list] = get_bedpe_header(args.bedpe)

	intersect = run_pair_to_bed(args.bedpe,args.regions,match_type="both").to_dataframe(header=None, comment="#")
	intersect = prep_intersection(intersect,bedpe_header_list)
	artefact_calls = get_artefact_calls(intersect)

	if args.intermediate:
			intersect.to_csv("intermediate.tsv",sep="\t",header=True, index=False)
	cdna_filter = pybedtools.BedTool(args.bedpe)
	cdna_filter = bedtool_to_df(cdna_filter,bedpe_header_list)
	cdna_filter = cdna_filter.merge(artefact_calls, on="ID", how="left")
		
	with open(args.outbedpe, "w") as fw:
		fw.write("".join(meta_header))
		
	cdna_filter.to_csv(args.outbedpe, header=True, index=False, sep="\t", mode="a")
	artefact_calls.to_csv(args.outfile, sep="\t", header=False, index=False)

if __name__ == "__main__":
	main()
