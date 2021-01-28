#!/opt/conda/envs/multiqc/bin/python3.9

__author__       = "Anne Marie Noronha"
__contributor__  = ""
__email__        = "noronhaa@mskcc.org"
__version__      = "0.0.2"

import json, os, yaml
import pandas as pd
from pandas.io.json import json_normalize
import sys, argparse
if sys.version_info[0] < 3: 
    from StringIO import StringIO
else:
    from io import StringIO

upperQuantile=.999
lowerQuantile=.001

#chromosomes = [str(i) for i in range(1,23)] + ["X","Y","MT"]
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
		result["_".join(["MQ",i])] = graphParser(files[i],i,"MQ","MappingQuality","Fraction") # linegraph if multiple samples
		result["_".join(["GC",i])] = graphParser(files[i],i,"GC","GCcontent","fractionOfReads",ignoreSamples=["Target","Reference"])
		result["_".join(["CO",i])] = graphParser(files[i],i,"CO","Coverage","Count") 
		result["_".join(["IS",i])] = graphParser(files[i],i,"IS","InsertSize","Count", addIDs=["Layout"])
		result["_".join(["IC",i])] = graphParser(files[i],i,"IC",["Homopolymer","InDel"],"Count") #bargraph
		result["_".join(["IC",i])] = defaultCPswitch(result["_".join(["IC",i])])
		for h in getHeaders(result["_".join(["IC",i])]):
			print("header: " + str(h))
			result["_".join(["IC",i])] = addColumnFormat(result["_".join(["IC",i])],h)
		result["_".join(["OT",i])] = graphParser(files[i],i,"OT","Extension","OnTarget") #linegraph
		result["_".join(["OT",i])] = addXCategory(result["_".join(["OT",i])])
		result["_".join(["CM",i])] = graphParser(files[i],i,"CM","Chrom","ObsExpRatio") #scatter
		result["_".join(["CM",i])] = addYPlotLine(addXCategory(result["_".join(["CM",i])]))
		result["_".join(["ME",i])] = MEParser(files[i],i) # generalstats

		printJSON("_".join(["IZ",i]),addTextLabels(result["_".join(["IZ",i])],i,"IZ","InDel Size","InDel Size by InDel type (INS or DEL)"))
		printJSON("_".join(["BQ",i]),addTextLabels(result["_".join(["BQ",i])],i,"BQ","Base Quality", "Base Quality by position in each read: Read1 and Read2"))
		printJSON("_".join(["MQ",i]),addTextLabels(result["_".join(["MQ",i])],i,"MQ","Mapping quality","Mapping quality by fraction of reads"))
		printJSON("_".join(["GC",i]),addTextLabels(result["_".join(["GC",i])],i,"GC","GC-Content","GC-Content by read count"))
		printJSON("_".join(["CO",i]),addTextLabels(result["_".join(["CO",i])],i,"CO","Coverage","Coverage by base in the reference"))
		printJSON("_".join(["IS",i]),addTextLabels(result["_".join(["IS",i])],i,"IS","Insert Size","Insert Size distribution"))
		printJSON("_".join(["IC",i]),addTextLabels(result["_".join(["IC",i])],i,"IC","Homopolymer distribution","Homopolymer distribution by InDel type (INS or DEL)","table"))
		printJSON("_".join(["OT",i]),addTextLabels(result["_".join(["OT",i])],i,"OT","On-Target rate","Rate of reads mapping within the target regions when each region is extended by a given extension length (bp) in either direction"))
		printJSON("_".join(["CM",i]),addTextLabels(result["_".join(["CM",i])],i,"CM","Chromosome Mapping","Observed vs expected ratio of mapping per chromosome, with y=1 plotted as a red line. Unplaced contigs ignored"))		
		printJSON("_".join(["ME",i]),result["_".join(["ME",i])])


def printJSON(label,jsonObj):
	#fname = label + '_mqc.json'
	fname = label + '_mqc.yaml'
	def invalidResult():
		print("No data for label " + label)
		if os.path.exists(fname):
			os.remove(fname)
	if jsonObj:
		if len(jsonObj) > 0:
			with open(fname, 'w') as outfile:
				#json.dump(jsonObj, outfile, indent=0)
				yaml.dump(jsonObj, outfile, indent=0)
			outfile.close()
		else: invalidResult()
	else: invalidResult()

def addPConfig(j,k,v):
	try:
		j["pconfig"][k] = v
	except:
		print("invalid object for adding pconfig")
	finally: return j

def addHeaderConfig(j,col,k,v):
	try:
		if "headers" not in j: 
			j["headers"] = dict()
		if col not in j["headers"]:
			j["headers"][col] = dict()
		j["headers"][col][k] = v
	except:
		print("invalid object for adding pconfig")
	finally: return j

def showCPswitch(jObj, showSwitch=False): #turns on the counts/percentages switch. the data should hold the counts data, percentages are calculated when generating the html. works with bargraphs, not linegraph
	return addPConfig(jObj,"cpswitch",showSwitch)

def toggle_xDecimals(jObj, xDecimal = False):
	return addPConfig(jObj,"xDecimals",xDecimal)

def toggle_yDecimals(jObj, yDecimal = False):
	return addPConfig(jObj,"yDecimals",yDecimal)

def addColumnFormat(jObj,column,format="{:,.0f}"):
	return addHeaderConfig(jObj,column,"format",format)

def defaultCPswitch(jObj,countByDefault=False): #sets the default positing of the counts/percentages switch. countByDefault=False makes the page load that section with Percentages view on. works with bargraphs, not linegraph
	return addPConfig(jObj,"cpswitch_c_active",countByDefault)

def addXPlotLine(jObj,value=1,color='#FF0000'):  # adds line at x = value
	return addPConfig(jObj,"xPlotLines",[{'color':color,'value':value}])
		
def addYPlotLine(jObj,value=1,color='#FF0000'): # adds line at y = value
	return addPConfig(jObj,"yPlotLines",[{'color':color,'value':value}])
	
def addXLog(jObj): #Use log10 x axis
	return addPConfig(jObj,"xLog",True)

def addYLog(jObj): #Use log10 y axis
	return addPConfig(jObj,"yLog",True)

def addXCategory(jObj): #Set to True to use x values as categories instead of numbers. 
	return addPConfig(jObj,"categories",True)

def getHeaders(jObj):
	retSet = set()
	try:
		retSet = set([ item for sub in [ list(jObj["data"][k].keys()) for k in jObj["data"] ] for item in sub ] )
	finally: return retSet

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
		ret_json["pconfig"]["title"] = title + " Read-group-" + RG
		ret_json["pconfig"]["id"] = sectionid + "_plot"
		ret_json["plot_type"] = plot_type
		ret_json["section_name"] = title + " read-group-" + RG
		ret_json["parent_id"] = "Alfred_QC_stats"
		ret_json["parent_name"] = 'Alfred QC'
		ret_json["parent_description"] = "Stats produced by <a href=\"https://www.gear-genomics.com/docs/alfred/\">AlfredQC</a> "
	except: 
		
		print("Something went wrong")
	finally: 
		return ret_json


def readFiles(listOfFiles,prefix):
	df = pd.DataFrame()
	cmd = "zgrep ^{} {}".format(prefix, " ".join(listOfFiles)) + " | cut -f 2- | sed '1!{/^Sample\\t/d;}'"
	filterFile = os.popen(cmd).read()
	try:
		df = pd.read_table(StringIO(filterFile), sep="\t", header=0)
	except:
		return None
	else:
		return df


def graphParser(listOfFiles,RG,prefix,x,y,addIDs=None,ignoreSamples=[],applyQuantiles=True):
	if len(listOfFiles) < 1: return None 
	ret_json = dict()
	xlab = x 
	ylab = y
	quantileCol = x	
	try:
		df = readFiles(listOfFiles,prefix)
		df = df[~df["Sample"].isin(ignoreSamples) ]
	except:
		print("something went wrong while filtering sample column (prefix {})".format(prefix))
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
		#ret_json[libraryID][str(row[xlab])] = row[ylab]
		ret_json[libraryID][row[xlab]] = row[ylab]
	j = list(ret_json.keys())
	for i in j:
		if len(ret_json[i]) == 0 or set(ret_json[i].values()) == {0}:
			ret_json.pop(i)
	if (len(ret_json) > 0):
		return {"pconfig":{"xlab":xlab,"ylab":ylab},"data":ret_json}
	else: return None

def MEParser(listOfFiles,RG,keepCols=["DuplicateFraction"]):
	if len(listOfFiles) < 1: return None 
	ret_json = dict()
	try:
		df = readFiles(listOfFiles,"ME")
		#df = df.filter(items=["Sample"].extend(keepCols))
		colPercents = [i for i in list(df) if "Fraction" in i or "Rate" in i ]
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
		#ret_json[libraryID] = { i:row[i]*100 if i in colPercents else row[i] for i in set(list(df)) - {"Sample","Library"} } 
		ret_json[libraryID] = { i:row[i]*100 if i in colPercents else row[i] for i in keepCols }
	#return {"plot_type":"generalstats","pconfig":[{i:{"suffix":"%" if i in colPercents else "" , "description": i + " extracted from Alfred QC", "hidden":not i in keepCols}} for i in list(df) if i not in ["Sample","Library"]],"data":ret_json}
	return {"plot_type":"generalstats","pconfig":[{i:{"suffix":"%" if i in colPercents else "" , "description": i + " extracted from Alfred QC", "hidden":False}} for i in list(df) if i in keepCols ],"data":ret_json}


if __name__ == "__main__":
	main()
