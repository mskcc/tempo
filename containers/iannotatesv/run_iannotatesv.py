import argparse, sys, os, subprocess
import numpy as np
import pandas as pd

def usage():
	parser = argparse.ArgumentParser()
	parser.add_argument('--bedpe', help = 'bedpe file', required = True)
	parser.add_argument('--genome', default = "hg19", help = 'hg19/hg38/hg18')

	return parser.parse_args()

def main():
	args = usage()

	meta=""
	with open(args.bedpe, 'r') as f:
		main_data = False
		while main_data == False:
			x = f.readline()
			if x.startswith("##"):
				meta += x
			else: 
				header=x
				main_data = True
		try:
			data = pd.read_csv(f, header=None, sep="\t" )
			data.columns = header.strip().split("\t")
		except:
			data = pd.DataFrame(columns = header.strip().split("\t"))

	iAnnotate_input = data["#CHROM_A|START_A|STRAND_A|CHROM_B|START_B|STRAND_B|ID".split("|")]
	iAnnotate_input.columns = "chr1|pos1|str1|chr2|pos2|str2|ID".split("|")

	iAnnotate_input = iAnnotate_input.replace("-", 1).replace("+", 0)

	iAnnotate_input.to_csv("input.txt",sep="\t",header=True,index=False)

	if not os.path.isdir("outputdir"):
		os.mkdir("outputdir")
	run_args = "python /usr/bin/iAnnotateSV/iAnnotateSV/iAnnotateSV.py -i input.txt -ofp outputfilePrefix -o outputdir -r {} -d 3000""".format(args.genome).split(" ")
	p = subprocess.Popen(run_args)
	p.communicate()

	iAnnotate_output = pd.read_csv("outputdir/outputfilePrefix_Annotated.txt", sep="\t",header=0)
	iAnnotate_output = pd.merge(iAnnotate_input,iAnnotate_output, on="chr1|pos1|str1|chr2|pos2|str2".split("|"), how="inner")

	annot_data = iAnnotate_output.drop("chr1|pos1|str1|chr2|pos2|str2".split("|"), axis=1)
	annot_data = annot_data.replace(np.nan, ".")

	final_data = pd.merge(data, annot_data, on="ID",how="left")

	outputbed = os.path.basename(args.bedpe)[:-6] if args.bedpe.endswith(".bedpe") else os.path.basename(args.bedpe)
	outputbed = os.path.join(os.path.dirname(args.bedpe),outputbed + ".iannotate.bedpe")
	with open(outputbed, "w") as fw:
		fw.write(meta)
	final_data.to_csv(outputbed,sep="\t",header=True,index=False, mode='a')

if __name__ == "__main__":
	main()
