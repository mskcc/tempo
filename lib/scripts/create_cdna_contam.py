#!/usr/bin/env python
import pandas as pd
import argparse

if __name__ == '__main__':

    parser = argparse.ArgumentParser(add_help=True, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--input_mafs', required=True, nargs='+', help='Need list of mafs, like *svs.pass.vep.maf')
    parser.add_argument('--project_prefix', required=True, help='Project name, ie Proj_000001_A')
    args = parser.parse_args()

    counter = 0
    for infile in args.input_mafs:
        print("current is fileo: " + infile)
        df = pd.read_csv(infile, sep='\t', comment='#')
        if 'vcf_id' in df.columns and 'Consequence' in df.columns and 'CONSENSUS' in df.columns:
            df['vartype'] = df['vcf_id'].apply(lambda x: x[0:3])
            deldf = df.loc[df['vartype']=='DEL'] #Grab only deletions
            deldf = deldf[deldf['Consequence'].str.contains("splice", case=False)] #Grab splice only
            deldf = deldf[pd.notnull(deldf['CONSENSUS'])] #Grab PRECISE
            if not deldf.empty:
                countdf = deldf.groupby(['Tumor_Sample_Barcode','Hugo_Symbol','vcf_id'])['vcf_id'].count().reset_index(name="count")
                countdf = countdf.loc[countdf['count'] >= 2] #make sure each variant_id has two instances
                if not countdf.empty:
                    cdnadf = countdf.groupby(['Tumor_Sample_Barcode', 'Hugo_Symbol'])['Hugo_Symbol'].count().reset_index(name="count")
                    if counter == 0:
                        counter+=1
                        cdnadf.to_csv('%s_cdna_contamination.txt' % args.project_prefix, mode='w', index=False, sep='\t')
                    else:
                        cdnadf.to_csv('%s_cdna_contamination.txt' % args.project_prefix, mode='a', header=False, index=False, sep='\t')
    if counter == 0:
        with open('%s_cdna_contamination.txt' % args.project_prefix, 'w') as tempfile:
            tempfile.write('No detected cDNA contamination for %s.\n' % args.project_prefix)
