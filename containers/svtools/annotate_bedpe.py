import argparse
import pandas as pd
from utils import *

keepCols = "ID|RE_gene|P_gain_phen|P_gain_hpo|P_gain_source|P_gain_coord|P_loss_phen|"\
"P_loss_hpo|P_loss_source|P_loss_coord|P_ins_phen|P_ins_hpo|P_ins_source|P_ins_coord|"\
"B_gain_source|B_gain_coord|B_gain_AFmax|B_loss_source|"\
"B_loss_coord|B_loss_AFmax|B_ins_source|B_ins_coord|B_ins_AFmax|B_inv_source|B_inv_coord|"\
"B_inv_AFmax|TAD_coordinate|ENCODE_experiment|GC_content_left|GC_content_right|"\
"Repeat_coord_left|Repeat_type_left|Repeat_coord_right|Repeat_type_right|Gap_left|"\
"Gap_right|SegDup_left|SegDup_right|ENCODE_blacklist_left|ENCODE_blacklist_characteristics_left|"\
"ENCODE_blacklist_right|ENCODE_blacklist_characteristics_right|HI|TS|ExAC_delZ|ExAC_dupZ|"\
"ExAC_cnvZ|ExAC_synZ|ExAC_misZ|"\
"LOEUF_bin|GnomAD_pLI|ExAC_pLI|AnnotSV_ranking_score|AnnotSV_ranking_criteria|ACMG_class"

keepCols=keepCols.split("|")
print(keepCols)

def usage():
	parser = argparse.ArgumentParser()
	parser.add_argument('--bp1', help = 'annotsv with first breakpoints result file', required = True)
	parser.add_argument('--bp2', help = 'annotsv with second breakpoints result file', required = True)
	#parser.add_argument('--annotsv', help = 'annotsv result file', required = True)
	parser.add_argument('--bedpe', help = 'bedpe file', required = True)
	parser.add_argument('--out',help = 'output file' , required = True)
	#parser.add_argument('--keepannotsv',help = 'annotsv columns to add to bedpe file' , required = True)
	return parser.parse_args()

def main():
	args = usage()
	bedpe = pybedtools.BedTool(args.bedpe)
	[meta_header, bedpe_header_list] = get_bedpe_header(args.bedpe)
	bedpe = bedtool_to_df(bedpe,bedpe_header_list)

	try:
		bp1 = pd.read_csv(args.bp1, sep="\t",header=0)
		bp2 = pd.read_csv(args.bp2, sep="\t",header=0)
	except Exception as e:
		print(e)
		bedpe = bedpe.reindex(columns = bedpe.columns.tolist() + ['GENE1','GENE2']) 
	else:
		bedpe = bedpe.merge(get_fusions(bp1,bp2), on="ID",how="left")

	with open(args.out, "w") as fw:
		fw.write("".join(meta_header))
		
	bedpe.to_csv(args.out, header=True, index=False, sep="\t", mode="a")

def get_fusions(bp1,bp2):
	cols = ["ID","Gene_name"]
	bp1 = bp1[(bp1.Annotation_mode == "split") & (bp1.Location2 == "CDS")][cols]
	bp2 = bp2[(bp2.Annotation_mode == "split") & (bp2.Location2 == "CDS")][cols]
	print(bp1.head())
	bp1 = bp1.groupby(['ID'], as_index = False).agg({'Gene_name': ';'.join}).reset_index(drop=True)
	bp2 = bp2.groupby(['ID'], as_index = False).agg({'Gene_name': ';'.join}).reset_index(drop=True)
	print(bp2.head())
	fusion = bp1.merge(bp2, how="outer",on="ID",suffixes=("","2")).rename(columns={"Gene_name":"GENE1","Gene_name2":"GENE2"})
	print(fusion.head())
	return fusion

if __name__ == "__main__":
	main()
