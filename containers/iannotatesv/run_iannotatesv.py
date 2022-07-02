import argparse, sys, os, subprocess
import numpy as np
import pandas as pd
from utils import *

def usage():
	parser = argparse.ArgumentParser()
	parser.add_argument('--bedpe', help = 'bedpe file', required = True)
	parser.add_argument('--genome', default = "hg19", help = 'hg19/hg38/hg18')

	return parser.parse_args()

def main():
	args = usage()

	print("Running iAnnotateSV on {}".format(os.path.basename(args.bedpe)))

	data = parse_svtools_bedpe_file(args.bedpe)
	iAnnotate_input = data["#CHROM_A|START_A|STRAND_A|CHROM_B|START_B|STRAND_B|ID".split("|")]
	iAnnotate_input.columns = "chr1|pos1|str1|chr2|pos2|str2|ID".split("|")

	iAnnotate_input = iAnnotate_input.replace("-", 1).replace("+", 0)

	n_chunks = int((999+iAnnotate_input.shape[0])/1000)
	iAnnotate_output = pd.DataFrame()

	for i in range(n_chunks):
		iAnnotate_input_chunk = iAnnotate_input.iloc[(i-1)*1000:i*1000]
		input_file = "input_" + str(i) + ".txt"
		output_dir = "output_" + str(i)
		iAnnotate_input_chunk.to_csv(input_file,sep="\t",header=True,index=False)
		if not os.path.isdir(output_dir): os.mkdir(output_dir)
		run_args = "python /usr/bin/iAnnotateSV/iAnnotateSV/iAnnotateSV.py -i {} -ofp outputfilePrefix -o {} -r {} -d 3000""".format(input_file, output_dir, args.genome).split(" ")
		p = subprocess.Popen(run_args)
		p.communicate()
		iAnnotate_output_chunk = pd.read_csv(os.path.join(output_dir,"outputfilePrefix_Annotated.txt"), sep="\t",header=0)
		iAnnotate_output_chunk = pd.merge(iAnnotate_input_chunk,iAnnotate_output_chunk, on="chr1|pos1|str1|chr2|pos2|str2".split("|"), how="inner")
		iAnnotate_output = pd.concat([iAnnotate_output,iAnnotate_output_chunk])

	annot_data = iAnnotate_output.drop("chr1|pos1|str1|chr2|pos2|str2".split("|"), axis=1)
	final_data = pd.merge(data, annot_data, on="ID",how="left")

	final_data = pd.merge(data, annot_data, on="ID",how="left")

	outputbed = os.path.basename(args.bedpe)[:-6] if args.bedpe.endswith(".bedpe") else os.path.basename(args.bedpe)
	outputbed = os.path.join(os.path.dirname(args.bedpe),outputbed + ".iannotate.bedpe")
	with open(outputbed, "w") as fw:
		fw.write(meta)
	final_data.to_csv(outputbed,sep="\t",header=True,index=False, mode='a')

if __name__ == "__main__":
	main()
