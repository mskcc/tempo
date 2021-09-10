#!/opt/conda/envs/multiqc/bin/python3.9
import json, sys

__author__       = "Anne Marie Noronha"
__contributor__  = ""
__email__        = "noronhaa@mskcc.org"
__version__      = "0.0.1"
	
def clean_fastp_json(inPath,outPath):
	with open(inPath, 'r') as data_file:
		data = json.load(data_file)
	keys = list(data.keys())
	for element in keys:
		if element in ['read1_after_filtering','read2_after_filtering'] :
			data.pop(element, None)
	keys = list(data['summary'].keys())
	for element in keys:
		if element in ['after_filtering'] :
			data['summary'].pop(element, None)
	with open(outPath, 'w') as data_file:
		json.dump(data, data_file)
	
if __name__ == "__main__":
	clean_fastp_json(sys.argv[1],sys.argv[2])
