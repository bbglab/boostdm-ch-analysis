# conda activate clonalh

# qmap submit chvariants_filter_vep.qmap

import pandas as pd
import os
from io import StringIO
import copy
import click
import subprocess
import gzip

@click.command()

@click.option('--input_maf_1',
              '-i2',
              required = True,
              help="path to gvcf file from uk biobank")

@click.option('--input_maf_2',
              '-i1',
              required = True,
              help="path to gvcf file from uk biobank")

@click.option('--output_folder',
              '-o',
              required = True,
              help="path to output folder")

def filter_chgenes(input_maf_1, input_maf_2, output_folder):

    ### 1. Open MAF file
    lines1 = ''.join([line for line in open(input_maf_1, 'rt') if not line.startswith("##")])
    ver_1 = pd.read_csv(StringIO(lines1), sep= '\t')
    lines2 = ''.join([line for line in open(input_maf_2, 'rt') if not line.startswith("##")])
    ver_2 = pd.read_csv(StringIO(lines2), sep= '\t')

    # Define function Eliminates "<*>" from ALT
    def eliminate_nonref(alt):
        if alt.split(",")[-1] == "<*>":
            wo_nonref = ",".join(alt.split(",")[:-1])
            return wo_nonref
        else:
            return alt

    ### 2. Obtain only alteration version 1
    # Split the information of each alteration
    split_ukbcol = ver_1[ver_1.columns[-1]].str.split(":",expand=True)
    split_ukbcol.columns = ['PL','DP','SP','ADF','ADR', 'AD']
    ver_1 = pd.concat([ver_1, split_ukbcol], axis=1)
    ver_1['AD_alt'] = ver_1['AD'].apply(lambda x: x.split(',')[1])
    ver_1['AD_ref'] = ver_1['AD'].apply(lambda x: x.split(',')[0])
    ver_1['ADF_alt'] = ver_1['ADF'].apply(lambda x: x.split(',')[1])
    ver_1['ADF_ref'] = ver_1['ADF'].apply(lambda x: x.split(',')[0])
    ver_1['ADR_alt'] = ver_1['ADR'].apply(lambda x: x.split(',')[1])
    ver_1['ADR_ref'] = ver_1['ADR'].apply(lambda x: x.split(',')[0])
    # Eliminate "<NON_REF>"
    ver_1["ALT"] = ver_1["ALT"].apply(lambda x: eliminate_nonref(x))
    ver_1_alt = ver_1[ver_1['ALT']!='']

    ### 2. Obtain only alteration version 2
    split_ukbcol = ver_2[ver_2.columns[-1]].str.split(":",expand=True)
    split_ukbcol.columns = ['PL','DP','SP','ADF','ADR', 'AD']
    ver_2 = pd.concat([ver_2, split_ukbcol], axis=1)
    ver_2['AD_alt'] = ver_2['AD'].apply(lambda x: x.split(',')[1])
    ver_2['AD_ref'] = ver_2['AD'].apply(lambda x: x.split(',')[0])
    ver_2['ADF_alt'] = ver_2['ADF'].apply(lambda x: x.split(',')[1])
    ver_2['ADF_ref'] = ver_2['ADF'].apply(lambda x: x.split(',')[0])
    ver_2['ADR_alt'] = ver_2['ADR'].apply(lambda x: x.split(',')[1])
    ver_2['ADR_ref'] = ver_2['ADR'].apply(lambda x: x.split(',')[0])
    # Eliminate "<NON_REF>"
    ver_2["ALT"] = ver_2["ALT"].apply(lambda x: eliminate_nonref(x))
    ver_2_alt = ver_2[ver_2['ALT']!='']

    ### Merge files
    final = pd.merge(ver_2, ver_1, left_index=True, right_index=True)
    #Obtain only mutation positions
    final = final.iloc[sorted(list(set(ver_1_alt.index.tolist()+ver_2_alt.index.tolist()))),:]
    # Merge reads
    def alt_merge(var):
        if (var['ALT_x'] == var['ALT_y']) :
            return var['ALT_x']
        elif (var['ALT_x'] == '') and (var['ALT_y']!=''):
            return var['ALT_y']
        elif (var['ALT_x'] != '') and (var['ALT_y']==''):
            return var['ALT_x']
        else:
            return var['ALT_x']+','+var['ALT_y']
    # Add type of variant

    if len(final)>0:
        final['ALT'] = final.apply(lambda x: alt_merge(x), axis=1)
        final['DP'] = final['DP_x'].astype(int)+final['DP_y'].astype(int)
        final['AD_alt'] = final['AD_alt_x'].astype(int)+final['AD_alt_y'].astype(int)
        final['AD_ref'] = final['AD_ref_x'].astype(int)+final['AD_ref_y'].astype(int)
        final['ADF_alt'] = final['ADF_alt_x'].astype(int)+final['ADF_alt_y'].astype(int)
        final['ADF_ref'] = final['ADF_ref_x'].astype(int)+final['ADF_ref_y'].astype(int)
        final['ADR_alt'] = final['ADR_alt_x'].astype(int)+final['ADR_alt_y'].astype(int)
        final['ADR_ref'] = final['ADR_ref_x'].astype(int)+final['ADR_ref_y'].astype(int)
        final['VAF_alt'] = final['AD_alt'].astype(int)/final['DP'].astype(int)

    	# Single allel
        final_single = final[~final['ALT'].str.contains(',')][['#CHROM_x', 'POS_x', 'ID_x', 'REF_x',\
                                                       'ALT', 'QUAL_x', 'FILTER_x', 'DP', 'AD_alt',\
       	'AD_ref', 'ADF_alt', 'ADF_ref','ADR_alt', 'ADR_ref', 'VAF_alt']]

    	# Multiple allel
        final_multi = final[final['ALT'].str.contains(',')]
        for n in final_multi.index.tolist():
            row = final.loc[n]
            alele_A = {'#CHROM_x':row['#CHROM_x'] , 'POS_x':row['POS_x'], 'ID_x':row['ID_x'], 'REF_x':row['REF_x'],\
                   	'ALT':row['ALT_x'], 'QUAL_x':row['QUAL_x'], 'FILTER_x':row['FILTER_x'],\
                   	'DP':row['DP'],\
                   	'AD_alt':row['AD_alt_x'],\
                   	'AD_ref':row['AD_ref'],\
                   	'ADF_alt':row['ADF_alt_x'],\
                   	'ADF_ref':row['ADF_ref'],\
                   	'ADR_alt':row['ADR_alt_x'],\
                   	'ADR_ref':row['ADR_ref'],\
                   	'VAF_alt':int(row['AD_alt_x'])/int(row['DP'])}
            final_single = final_single.append(alele_A, ignore_index=True)

            alele_B = {'#CHROM_x':row['#CHROM_x'] , 'POS_x':row['POS_x'], 'ID_x':row['ID_x'], 'REF_x':row['REF_x'],\
                   	'ALT':row['ALT_y'], 'QUAL_x':row['QUAL_x'], 'FILTER_x':row['FILTER_y'],\
                   	'DP':row['DP'],\
                   	'AD_alt':row['AD_alt_y'],\
                   	'AD_ref':row['AD_ref'],\
                   	'ADF_alt':row['ADF_alt_y'],\
                   	'ADF_ref':row['ADF_ref'],\
                   	'ADR_alt':row['ADR_alt_y'],\
                   	'ADR_ref':row['ADR_ref'],\
                   	'VAF_alt':int(row['AD_alt_y'])/int(row['DP'])}
            final_single = final_single.append(alele_B, ignore_index=True)

    else:
        final['ALT'] = []
        final['ALT'] = []
        final['DP'] = []
        final['AD_alt'] = [] 
        final['AD_ref'] = []
        final['ADF_alt'] = []
        final['ADF_ref'] = []
        final['ADR_alt'] = []
        final['ADR_ref'] = []
        final['VAF_alt'] = []
        final_single = final[['#CHROM_x', 'POS_x', 'ID_x', 'REF_x',\
                             'ALT', 'QUAL_x', 'FILTER_x', 'DP', 'AD_alt',\
       			'AD_ref', 'ADF_alt', 'ADF_ref','ADR_alt', 'ADR_ref', 'VAF_alt']]

    final_single['VEP'] = ''
    final_single = final_single[['#CHROM_x', 'POS_x', 'ID_x', 'REF_x', 'ALT', 'QUAL_x', 'FILTER_x',\
                             'VEP', 'AD_alt', 'DP', 'VAF_alt', 'ADF_alt', 'ADR_alt']]

    ### Save file
    # Get file name without extension form initial path and create temp output folder+file
    file = input_maf_1.split("/")[-1].split("_2.")[0]
    temp_file = output_folder + file + ".maf.gz"
    # Save
    final_single.to_csv(temp_file, sep="\t", index=False, compression='gzip')

if __name__ == "__main__":
    filter_chgenes()
