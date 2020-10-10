#!/opt/conda/envs/multiqc/bin/python3.8

__author__       = "Anne Marie Noronha"
__contributor__  = ""
__email__        = "noronhaa@mskcc.org"
__version__      = "0.0.1"

import json, os
import pandas as pd
from pandas.io.json import json_normalize
import sys, argparse
if sys.version_info[0] < 3: 
    from StringIO import StringIO
else:
    from io import StringIO

alfred_prefix={"Coverage":"CO",
	       "Mapping Quality":"MQ",
	       "Read Length":"RL",
	       "Base Qualities":"BQ",
	       "Base Content":"BC",
	       "Insert Size":"IS",
	       "InDel Context":"IC",
	       "Chromosome GC-content":"CG",
	       "Avg. target coverage":"TC",
	       "On target rate":"OT",
	       "Alignment summary metrics":"ME",
	       "Chromosome mapping statistics":"CM",
	       "InDel size":"IZ",
	       "GC-content":"GC"
}

graphtype={ i:"linegraph" for i in alfred_prefix}
graphtype["Alignment summary metrics"] = "table"

upperQuantile=.999
lowerQuantile=.001

chromosomes = [str(i) for i in range(1,23)] + ["X","Y","MT"]
chromosomes = [str(i) for i in range(1,23)] + ["X","Y"]

def usage():
	parser = argparse.ArgumentParser(description='Search Lims according to Request ID or Sample ID')
	parser.add_argument('--alfredfiles', dest='alfredFileList', metavar='XXXXX.tsv.gz', type=str, nargs='+', required=False,
	                    help='specify as many request IDs as desired separated by spaces')
	parser.add_argument('--multiple-samples',action="store_true",
						help='specify whether to produce single sample view or multi-sample view')
	args = parser.parse_args()
	return args

def main():
	args = usage()
	rgawarefiles  = [ x for x in args.alfredFileList if "per_readgroup" in x ]
	rgignorefiles = [ x for x in args.alfredFileList if "per_readgroup" not in x ]
	files = {"aware":rgawarefiles, "ignore":rgignorefiles}

	result = {}

	for i in ["aware","ignore"]:
		result["_".join(["IZ",i])] = graphParser(files[i],i,"IZ","Size","Count",addIDs=["InDel"])
		result["_".join(["BQ",i])] = graphParser(files[i],i,"BQ","Position","BaseQual",addIDs=["Read"])
		result["_".join(["MQ",i])] = graphParser(files[i],i,"MQ","MappingQuality","Count") # linegraph if multiple samples
		result["_".join(["GC",i])] = graphParser(files[i],i,"GC","GCcontent","fractionOfReads",ignoreSamples=["Target","Reference"])
		result["_".join(["CO",i])] = graphParser(files[i],i,"CO","Coverage","Count") 
		result["_".join(["IS",i])] = graphParser(files[i],i,"IS","InsertSize","Count", addIDs=["Layout"])
		result["_".join(["IC",i])] = graphParser(files[i],i,"IC",["Homopolymer","InDel"],"Count") #bargraph
		result["_".join(["IC",i])] = defaultCPswitch(result["_".join(["IC",i])])
		result["_".join(["OT",i])] = graphParser(files[i],i,"OT","Extension","OnTarget") #linegraph
		result["_".join(["OT",i])] = addXCategory(result["_".join(["OT",i])])
		result["_".join(["CM",i])] = graphParser(files[i],i,"CM","Chrom","ObsExpRatio") #scatter
		result["_".join(["CM",i])] = addYPlotLine(addXCategory(result["_".join(["CM",i])]))
		result["_".join(["ME",i])] = MEParser(files[i],i) # generalstats

		printJSON("_".join(["IZ",i]),addTextLabels(result["_".join(["IZ",i])],i,"IZ","InDel Size","InDel Size by InDel type (INS or DEL)"))
		printJSON("_".join(["BQ",i]),addTextLabels(result["_".join(["BQ",i])],i,"BQ","Base Quality", "Base Quality by position in each read: Read1 and Read2"))
		printJSON("_".join(["MQ",i]),addTextLabels(result["_".join(["MQ",i])],i,"MQ","Mapping quality","Mapping quality by read count"))
		printJSON("_".join(["GC",i]),addTextLabels(result["_".join(["GC",i])],i,"GC","GC-Content","GC-Content by read count"))
		printJSON("_".join(["CO",i]),addTextLabels(result["_".join(["CO",i])],i,"CO","Coverage","Coverage by base in the reference"))
		printJSON("_".join(["IS",i]),addTextLabels(result["_".join(["IS",i])],i,"IS","Insert Size","Insert Size distribution"))
		printJSON("_".join(["IC",i]),addTextLabels(result["_".join(["IC",i])],i,"IC","Homopolymer distribution","Homopolymer distribution by InDel type (INS or DEL)","bargraph"))
		printJSON("_".join(["OT",i]),addTextLabels(result["_".join(["OT",i])],i,"OT","On-Target rate","On-Target rate"))
		printJSON("_".join(["CM",i]),addTextLabels(result["_".join(["CM",i])],i,"CM","Chromosome Mapping","Observed vs expected ratio of mapping per chromosome, with y=1 plotted as a red line. Unplaced contigs ignored"))		
		printJSON("_".join(["ME",i]),result["_".join(["ME",i])])		


def printJSON(label,jsonObj):
	fname = label + '_mqc.json'
	def invalidResult():
		print("No data for label " + label)
		if os.path.exists(fname):
			os.remove(fname)
	if jsonObj:
		if len(jsonObj) > 0:
			with open(fname, 'w') as outfile:
				json.dump(jsonObj, outfile, indent=0)
			outfile.close()
		else: invalidResult()
	else: invalidResult()

def addPConfig(j,k,v):
	try:
		j["pconfig"][k] = v
	except:
		print("invalid object for adding pconfig")
	finally: return j

def showCPswitch(jObj, showSwitch=False):
	return addPConfig(jObj,"cpswitch",showSwitch)

def defaultCPswitch(jObj,on_by_default=False):
	return addPConfig(jObj,"cpswitch_c_active",on_by_default)

def addXPlotLine(jObj,value=1,color='#FF0000'):
	return addPConfig(jObj,"xPlotLines",[{'color':color,'value':value}])
		
def addYPlotLine(jObj,value=1,color='#FF0000'):
	return addPConfig(jObj,"yPlotLines",[{'color':color,'value':value}])
	
def addXLog(jObj):
	return addPConfig(jObj,"xLog",True)

def addYLog(jObj):
	return addPConfig(jObj,"yLog",True)

def addXCategory(jObj):
	return addPConfig(jObj,"categories",True)
		
def addTextLabels(base, RG, prefix, title,description,plot_type="linegraph"):
	ret_json = base
	try:
		if "pconfig" not in ret_json:
			ret_json["pconfig"] = dict()
		sectionid = "_".join(["Alfred", prefix.strip(),RG])
		ret_json["description"] = description + ", extracted from <a href=\"https://www.gear-genomics.com/docs/alfred/\">AlfredQC</a> result."
		if RG == "aware":
			ret_json["description"] = ret_json["description"] + " Showing results per read-group."
		ret_json["id"] = sectionid
		ret_json["pconfig"]["title"] = "Alfred " + title + " Read-group-" + RG
		ret_json["pconfig"]["id"] = sectionid + "_plot"
		ret_json["plot_type"] = plot_type
		ret_json["section_name"] = "Alfred " + title + " RG-" + RG
		ret_json["parent_id"] = "Alfred_QC_stats"
		ret_json["parent_name"] = 'Alfred QC'
		ret_json["parent_description"] = "Stats produced by <a href=\"https://www.gear-genomics.com/docs/alfred/\">AlfredQC</a> "
	except: 
		
		print("Something went wrong")
	finally: 
		#print(ret_json)
		return ret_json


def readFiles(listOfFiles,prefix):
	df = pd.DataFrame()
	cmd = "zgrep ^{} {}".format(prefix, " ".join(listOfFiles)) + " | cut -f 2- | sed '1!{/^Sample\\t/d;}'"
	filterFile = os.popen(cmd).read()
	df = pd.read_table(StringIO(filterFile), sep="\t", header=0)
	#for i in listOfFiles:
	#	cmd = "zgrep ^{} {} | cut -f 2- ".format(prefix, i)
	#	filterFile = os.popen(cmd).read()
	#	df = pd.concat([df,pd.read_table(StringIO(filterFile), sep="\t", comment="#", header=0)], ignore_index=True)
	return df


def graphParser(listOfFiles,RG,prefix,x,y,addIDs=None,ignoreSamples=[],applyQuantiles=True):
	if len(listOfFiles) < 1: return None 
	ret_json = dict()
	xlab = x 
	ylab = y
	quantileCol = x	
	df = readFiles(listOfFiles,prefix)
	try:
		df = df[~df["Sample"].isin(ignoreSamples) ]
	except:
		print("something went wrong while filtering sample column")
		return None
	if isinstance(x, list):
		df["__".join(x)] = ["@".join(i) for i in zip(*[df[x[i]].map(str) for i in range(len(x))])]
		xlab = "__".join(x)
		quantileCol = x[0]
	if "Quantile" in list(df) and applyQuantiles:
		df = df[(df["Quantile"] >= lowerQuantile) & (df["Quantile"] <= upperQuantile)]
	if "Chrom" in list(df):
		df = df[df["Chrom"].isin(chromosomes)]
	for index, row in df.iterrows():
		idLabels = [row["Sample"]]
		if RG == "aware": idLabels = [row["Sample"], row["Library"]]
		if addIDs: 
			idLabels.extend([row[i] for i in addIDs])
		libraryID= "@".join(idLabels)
		if libraryID not in ret_json:
			ret_json[libraryID] = dict()
		ret_json[libraryID][str(row[xlab])] = row[ylab]
	j = list(ret_json.keys())
	for i in j:
		if len(ret_json[i]) == 0 or set(ret_json[i].values()) == {0}:
			ret_json.pop(i)
	if (len(ret_json) > 0):
		return {"pconfig":{"xlab":xlab,"ylab":ylab,"xDecimals":False},"data":ret_json}
	else: return None

def CMParser(listOfFiles,RG,prefix,x="Chrom",y="ObsExpRatio",addIDs=None,ignoreSamples=[],applyQuantiles=True):
	if len(listOfFiles) < 1: return None 
	ret_json = dict()
	xlab = x 
	ylab = y
	df = readFiles(listOfFiles,prefix)
	try:
		df = df[~df["Sample"].isin(ignoreSamples) ]
	except:
		print("something went wrong while filtering sample column")
		return None
	try:
		df = df[df["Chrom"].isin(chromosomes)]
	except:
		print("This file does not have Chrom info")
		return None
	if isinstance(x, list):
		df["__".join(x)] = ["@".join(i) for i in zip(*[df[x[i]].map(str) for i in range(len(x))])]
		xlab = "__".join(x)
	for index, row in df.iterrows():
		idLabels = [row["Sample"]]
		if RG == "aware": idLabels = [row["Sample"], row["Library"]]
		if addIDs: 
			idLabels.extend([row[i] for i in addIDs])
		libraryID= "@".join(idLabels)
		if libraryID not in ret_json:
			ret_json[libraryID] = dict()
		ret_json[libraryID][str(row[xlab])] = row[ylab]
	j = list(ret_json.keys())
	for i in j:
		if len(ret_json[i]) == 0 or set(ret_json[i].values()) == {0}:
			ret_json.pop(i)
	if (len(ret_json) > 0):
		return {"pconfig":{"xlab":xlab,"ylab":ylab,"xDecimals":False},"data":ret_json}
	else: return None

def MEParser(listOfFiles,RG,keepCols=["DuplicateFraction"]):
	ret_json = dict()
	try:
		df = readFiles(listOfFiles,"ME")
		#df = df.filter(items=["Sample"].extend(keepCols))
		df = df.rename(columns={i:i.replace("#","Num") for i in list(df) if "#" in i})
	except:
		print("unable to read ME lines")
	for index, row in df.iterrows():
		print(row["Sample"])
		print(row["Library"])
		if RG == "aware": 
			libraryID = "@".join([row["Sample"], row["Library"]])
		else:
			libraryID = row["Sample"]
		ret_json[libraryID] = { i:row[i] for i in list(df) if i not in ["Sample","Library"]} 
	return {"plot_type":"generalstats","pconfig":[{i:{"description": i + " extracted from Alfred QC", "hidden":not i in keepCols}} for i in list(df) if i not in ["Sample","Library"]],"data":ret_json}


if __name__ == "__main__":
	main()
