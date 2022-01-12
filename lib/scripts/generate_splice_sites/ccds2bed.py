import pandas as pd
import sys
import urllib, os


url="https://ftp.ncbi.nlm.nih.gov/pub/CCDS/archive/Hs37.3/CCDS.current.txt"
source_file=os.path.basename(url)

if not os.path.exists(source_file):
	print("downloading source ccds file")
	urllib.urlretrieve(url,source_file)

ccds_df = pd.read_csv(source_file,sep="\t",header=0)
ccds_df = ccds_df[~ccds_df.ccds_status.str.contains("Withdrawn") & ~ccds_df.match_type.str.contains("Partial")]
ccds_df["cds_locations"] = ccds_df.apply(lambda x: [(k,v) for (k,v) in enumerate(x["cds_locations"].split(','))] if x["cds_strand"] == "+" else [(k,v) for (k,v) in enumerate(x["cds_locations"].split(',')[::-1])][::-1], axis =1)
ccds_df = ccds_df.explode('cds_locations').reset_index(drop=True)
ccds_df["exon_number"] = ccds_df.cds_locations.apply(lambda x:x[0])
ccds_df["cds_locations"] = ccds_df.apply(lambda x: list(zip(["acceptor","donor"] if x["cds_strand"] == "+" else ["donor","acceptor"],x["cds_locations"][1].strip("[] ").split("-") )),axis=1)
ccds_df = ccds_df.explode('cds_locations').reset_index(drop=True)
ccds_df["splice_orientation"] = ccds_df.cds_locations.apply(lambda x: x[0])
ccds_df["cds_locations"] = ccds_df.cds_locations.apply(lambda x: int(x[1]))
ccds_df = ccds_df.astype({"cds_from":"int64", "cds_to":"int64"})
ccds_df = ccds_df[(ccds_df.cds_locations != ccds_df.cds_from) & (ccds_df.cds_locations != ccds_df.cds_to)]
ccds_df = ccds_df["#chromosome|gene|ccds_id|cds_strand|cds_locations|exon_number|splice_orientation".split("|")]
ccds_df["start"] = ccds_df.apply(lambda x: x["cds_locations"] - 3 if (x["cds_strand"] == "+" and x["splice_orientation"] == "acceptor") or (x["cds_strand"] == "-" and x["splice_orientation"] == "donor") else x["cds_locations"] - 8,axis=1)
ccds_df["end"]   = ccds_df.apply(lambda x: x["cds_locations"] + 8 if (x["cds_strand"] == "+" and x["splice_orientation"] == "acceptor") or (x["cds_strand"] == "-" and x["splice_orientation"] == "donor") else x["cds_locations"] + 3,axis=1)
ccds_df["name"]  = ccds_df.apply(lambda x: "|".join([i+":"+str(x[i]) for i in "exon_number|splice_orientation|gene|ccds_id".split("|")]),axis=1)
ccds_df["score"] = '.'

ccds_df= ccds_df["#chromosome|start|end|name|score|cds_strand".split("|")]

ccds_df.to_csv("splice_sites.bed",header=False, index=False, sep="\t")