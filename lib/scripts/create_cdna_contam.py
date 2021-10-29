#!/usr/bin/env python
import pandas as pd
import argparse, sys, os

__author__  = "DMP"
__contributor__  = "Anne Marie Noronha"
__email__   = "noronhaa@mskcc.org"
__version__ = "0.0.1"
#__status__  = "Dev"


def usage():
    parser = argparse.ArgumentParser(add_help=True, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    #parser.add_argument('--input_mafs', required=True, nargs='+', help='Need list of mafs, like *svs.pass.vep.maf')
    #parser.add_argument('--project_prefix', required=True, help='Project name, ie Proj_000001_A')
    parser.add_argument('--input_maf', required=True, type=str, help='single sample maf file')
    parser.add_argument('--prefix', default='sample', type=str, help='prefix of output file')
    return parser.parse_args()
    
def read_splice_deletions(in_maf):
    try:
        df = pd.read_csv(in_maf, sep='\t', comment='#')
        df['vartype'] = df['vcf_id'].apply(lambda x: x[0:3])
        deldf = df.loc[df['vartype']=='DEL'] #Grab only deletions
        deldf = deldf[deldf['Consequence'].str.contains("splice", case=False)] #Grab splice only
        deldf = deldf[pd.notnull(deldf['CONSENSUS'])] #Grab PRECISE
    except:
        raise IOError("Unable to read and filter input maf. Check columns vcf_id, Consequence, CONSENSUS.")
    else:
        return deldf

def count_cdna_variants(deldf):
    try:
        countdf = deldf.groupby(['Tumor_Sample_Barcode','Hugo_Symbol','vcf_id'])['vcf_id'].count().reset_index(name="count")
        countdf = countdf.loc[countdf['count'] >= 2] #make sure each variant_id has two instances
        cdnadf = countdf.groupby(['Tumor_Sample_Barcode', 'Hugo_Symbol'])['Hugo_Symbol'].count().reset_index(name="count")
    except:
        raise ValueError("Unable to count cDNA variants. Check columns Tumor_Sample_Barcode, Hugo_Symbol, vcf_id")
    else: 
        return cdnadf

def output_cdna_counts(cdnadf, prefix):
    out_path = '%s_cdna_contamination.txt' % prefix
    if os.path.exists(out_path):
        print("appending genes..")
        mode = 'a'
        header=False
    else:
        print("writing to new file..")
        mode = 'w'
        header=True
        
    cdnadf.to_csv(out_path, mode=mode, header=header, index=False, sep='\t')


def main():
    args = usage()
    try:
        rawdf = read_splice_deletions(args.input_maf)
        countdf = count_cdna_variants(rawdf)
        output_cdna_counts(countdf, args.prefix)
    except (ValueError,IOError) as e:
        print(e)


if __name__ == '__main__':
    main()

