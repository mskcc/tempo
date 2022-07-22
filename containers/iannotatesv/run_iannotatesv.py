import argparse, sys, os, subprocess
from collections import OrderedDict
from datetime import datetime
import multiprocessing as mp
import numpy as np
import pandas as pd
from utils import *

def usage():
	parser = argparse.ArgumentParser()
	parser.add_argument('--bedpe', help = 'bedpe file', required = True)
	parser.add_argument('--genome', default = "hg19", help = 'hg19/hg38/hg18')
	parser.add_argument('--threads', type=int, default = 4, help = 'number of threads for parallelization')

	return parser.parse_args()

def run_iannotate_cmd(input_file,output_dir,output_pre,genome="hg19"):
	print("Spawning in parallel: annotation of " + input_file)
	run_args = "python /usr/bin/iAnnotateSV/iAnnotateSV/iAnnotateSV.py -i {} -o {} -ofp {} -r {} -d 3000""".format(input_file, output_dir, output_pre, genome).split(" ")
	p = subprocess.Popen(run_args)
	p.communicate()

def read_iannotatesv_result(path):
	return pd.read_csv(path, sep="\t",header=0)

def main():
	args = usage()

	dt = datetime.now()
	print("[{}] Running iAnnotateSV on {}".format(dt,os.path.basename(args.bedpe)))

	meta, header_list, data = parse_svtools_bedpe_file(args.bedpe)
	iAnnotate_input = data["#CHROM_A|START_A|STRAND_A|CHROM_B|START_B|STRAND_B|ID".split("|")]
	iAnnotate_input = iAnnotate_input.replace("-", 1).replace("+", 0)
	key_col = OrderedDict(zip("chr1|pos1|str1|chr2|pos2|str2".split("|"), [str,int,int,str,int,int]))
	iAnnotate_input.columns = key_col.keys() + ["ID"]
	for k,v in key_col.items():
		iAnnotate_input[k] = iAnnotate_input[k].astype(v)

	k = int(500)
	n_chunks = int(((k-1)+iAnnotate_input.shape[0])/k)
	iAnnotate_output = pd.DataFrame()
	if not os.path.isdir("inputs"): os.mkdir("inputs")
	if not os.path.isdir("outputs"): os.mkdir("outputs")

	pool = mp.Pool(args.threads)
	for i in range(n_chunks):
		iAnnotate_input_chunk = iAnnotate_input.iloc[i*k:(i+1)*k]
		input_file = os.path.join("inputs","input_" + str(i) + ".txt")
		iAnnotate_input_chunk.to_csv(input_file,sep="\t",header=True,index=False)
		pool.apply_async(run_iannotate_cmd, args=(input_file, "outputs", "output_" + str(i), args.genome ))

	pool.close()
	pool.join()

	iAnnotate_output = pd.concat(map(read_iannotatesv_result, [os.path.join("outputs","output_" + str(i) + "_Annotated.txt") for i in range(n_chunks)]))
	for key,val in key_col.items():
		iAnnotate_output[key] = iAnnotate_output[key].astype(val)

	iAnnotate_output = pd.merge(iAnnotate_input,iAnnotate_output, on=key_col.keys(), how="inner")

	print("[{}] Merging annotation with bedpe".format(datetime.now()))
	annot_data = iAnnotate_output.drop(key_col.keys(), axis=1).replace(np.nan, ".")
	final_data = pd.merge(data, annot_data, on="ID",how="left")

	outputbed = os.path.basename(args.bedpe)[:-6] if args.bedpe.endswith(".bedpe") else os.path.basename(args.bedpe)
	outputbed = os.path.join(os.path.dirname(args.bedpe),outputbed + ".iannotate.bedpe")
	with open(outputbed, "w") as fw:
		fw.write(meta)
	final_data.to_csv(outputbed,sep="\t",header=True,index=False, mode='a')

	print("[{}] iAnnotateSV complete".format(datetime.now()))

if __name__ == "__main__":
	main()
