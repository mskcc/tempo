import argparse, sys, os
import numpy as np
import pandas as pd
import pybedtools
from utils import *

def get_artefact_calls(df):
	"""
	Find genes with at least 2 exon-exon junctions
	Input requires SVTYPE, EXON, Hugo_Symbol
	Input requires variant id, gene names, exon names, splice site type
	remove any variant that doesn't involve two exons
	remove any variant that doesn't have "donor" and "acceptor" 
	remove any variant that doesn't have "donor" and "acceptor" in the correct orientation
	Only genes with at least two splice variants will be considered. 
	"""
	df = df.groupby(["var_id","gene"]).filter(lambda x: len(set(x['exon'])) == 2 and len(set(x['ss_type'])) == 2)
	#df = df.groupby(["var_id","gene"]).filter(lambda x: "donor" == sorted(zip(x['exon_number'],x['ss_type']), key = lambda t:t[1], reverse= True if x['strand'] == "-" else False)[0][1] )
	df = df.groupby(["var_id","gene"]).filter(lambda x: is_orientation_correct(x))
	df = df.groupby(["gene"]).filter(lambda x: x.size() > 1)
	return df["var_id"].tolist()
	
	#df = df.drop_duplicates().groupby(["var_id","gene"], as_index = False).agg({'exon': lambda x: len(set(x)) })
	#artefact_genes = df[df.exon > 3].gene.tolist()
	#return df[df.gene.isin(artefact_genes)]["var_id"].tolist()

def is_orientation_correct(x):
	reverse_ = True if x.strand.tolist()[0] == "-" else False
	zipped = sorted(zip(x['exon_number'],x['ss_type']), key = lambda t:t[0], reverse= reverse_)
	if zipped[0][1] == "donor":
		return True
	else: return False

def usage():
	parser = argparse.ArgumentParser()
	parser.add_argument('--exon-junct', dest="regions", help = 'regions with exon junctions', required = True)
	parser.add_argument('--bedpe', help = 'bedpe file', required = True)
	parser.add_argument('--output',dest="outfile",help = 'output file' , required = True)
    
	return parser.parse_args()

def main():
	args = usage()

	[meta_header, bedpe_header_list] = get_bedpe_header(args.bedpe)
	svtype_col = bedpe_header_list.index("TYPE")

	intersect = run_pair_to_bed(args.bedpe,args.regions,match_type="both")
	exon_del_df = pd.DataFrame(columns=["var_id","exon","exon_number","gene","ss_type","strand"])
	for i in intersect:
		if i[svtype_col] != "DEL":
			continue
		print(i)
		data_dict = dict()
		data_dict["var_id"] = i[6]
		data_dict["exon"] = i[ (len(bedpe_header_list) - 1) + 4 ].split(":")[0]
		data_dict["exon_number"] = data_dict["exon"].split(".")[-1]
		data_dict["gene"] = ".".join(data_dict["exon"].split(".")[:-1])
		# acceptor or donor
		data_dict["ss_type"] = i[ (len(bedpe_header_list) - 1) + 4 ].split(":")[1] 
		data_dict["strand"] = i[ (len(bedpe_header_list) - 1) + 6 ]
		exon_del_df = pd.concat([ exon_del_df, pd.DataFrame([data_dict]) ])

	exon_del_df.to_csv("debug.tsv",sep="\t",header=True, index=False)

	artefact_calls = get_artefact_calls(exon_del_df.drop_duplicates(keep='first'))

	cdna_filter = pybedtools.BedTool(args.bedpe)
	cdna_filter = bedtool_to_df(cdna_filter,bedpe_header_list)
	cdna_filter = cdna_filter.apply(lambda x: update_filter(x,"cdna") if x["ID"] in artefact_calls else x, axis=1)
	
	with open(args.outfile, "w") as fw:
		fw.write("".join(meta_header))
		
	cdna_filter.to_csv(args.outfile, header=True, index=False, sep="\t", mode="a")

if __name__ == "__main__":
	main()
