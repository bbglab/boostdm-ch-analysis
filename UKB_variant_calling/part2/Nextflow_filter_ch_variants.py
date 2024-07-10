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

@click.option('--input_maf',
              '-i',
              required = True,
              help="path to gvcf file from uk biobank")
@click.option('--output_folder',
              '-o',
              required = True,
              help="path to output folder")

def filter_chgenes(input_maf, output_folder):
    """
    Starting from MAF files from UKB 50K (gVCF filterd by CH gene coordinates)
    filters and arranges the data: select real variants, split variants with 2+ ALTs,
    annotate the data, etc AND run VEP
    """

    ### 1. Open MAF file
    lines = ''.join([line for line in gzip.open(input_maf, 'rt') if not line.startswith("##")])
    maf = pd.read_csv(StringIO(lines), sep= '\t')

    ### 2. Rearrange info from last column
    # Split column and add column names
    split_ukbcol = maf[maf.columns[-1]].str.split(":",expand=True).iloc[:,0:4]
    split_ukbcol.columns = ["GT","AD","AF","DP"]
    #Select columns from original maf and merge with splitted columns
    maf_var_spl = pd.concat([maf.iloc[:,[0,1,2,3,4,5,6]],
                                  split_ukbcol,
                                  maf.iloc[:,[7,8]]], axis=1)

    ### 3. Eliminate variants with AD=None, DP=None, or DP=0
    maf_var_spl = maf_var_spl[~maf_var_spl['INFO'].str.contains('DP=0')]
    maf_var_spl = maf_var_spl[~((maf_var_spl['AD'].isnull()) |
                            (maf_var_spl['DP'].isnull()) |
                            (maf_var_spl['AD'].str.split(',', expand=True)[0].isnull()) |
                            (maf_var_spl['AD'].str.split(',', expand=True)[1].isnull()) |
                            (maf_var_spl['DP']=="0"))]

    ### 4. Split variants with >1 ALT & Calculate VAF

    ### 4A. Variants 1 ALT
    # Select variants with 1 ALT
    df1 = maf_var_spl[~maf_var_spl['ALT'].str.contains(',')]
    # Calculate VAF (& add number of ALT=1)
    df1['AD_alt'] = df1['AD'].str.split(',', expand=True)[1].astype(int)
    df1['VAF_alt'] = df1['AD'].str.split(',', expand=True)[1].astype(int) / df1['DP'].astype(int)
    df1['VAF_ref'] = df1['AD'].str.split(',', expand=True)[0].astype(int) / df1['DP'].astype(int)
    df1['ALT_num'] = 1

    ### 4B. Variants 2+ ALT
    # Select variants with 2+ ALT & transform df to list
    df2 = maf_var_spl[maf_var_spl['ALT'].str.contains(',')]
    df2_list = df2.values.tolist()
    if len(df2_list)>0:
        #Divide variants in 1 line per ALT & calculate VAF
        df2_newlist = []
        # Loop through all variants
        for row in df2_list:
            # Extract info of ALTs and VAF
            ALTs = row[4].split(',')
            ADs = row[8].split(',')
            DP = row[10]
            # Loop to create 1 line per ALT
            for i in range(0,len(ALTs)):
                newrow = copy.deepcopy(row)
                # Take ALT and substitute column
                newrow[4] = ALTs[i]
                # Take AD corresponding ALT in new column
                newrow.append([])
                newrow[13] = int(ADs[i+1])
                # Calculate VAF from corresponding ALT in new column
                newrow.append([])
                newrow[14] = int(ADs[i+1]) / int(DP)
                # Calculate VAF from REF in new column
                newrow.append([])
                newrow[15] = int(ADs[0]) / int(DP)
                # Annotate number of total ALT in new column
                newrow.append([])
                newrow[16] = len(ALTs)
                # Append variant to new list
                df2_newlist.append(newrow)
                df2_newdf = pd.DataFrame(df2_newlist)
                df2_newdf.columns = df1.columns
                # ### 4C. Concatenate variants 1ALT and 2+ALT
                maf_var_spl_1alt = pd.concat([df1, df2_newdf], ignore_index=True)
    else:
        maf_var_spl_1alt = df1

    ### 5. Annotate variant type
    # Define function
    def variant_type(var):
        if (len(var[3]) == len(var[4])) & (len(var[3]) == 1):
            return 'SNV'
        elif len(var[3]) != len(var[4]):
            return 'Indel'
        elif (len(var[3]) == len(var[4])) & (len(var[3]) > 1):
            return 'MNV'
        else:
            return 'Unknown'
    # Add type of variant
    maf_var_spl_1alt['var_type'] = maf_var_spl_1alt.apply(lambda x: variant_type(x), axis=1)


    ### 6. Reorder columns and create new column for VEP output
    maf_var_spl_1alt.iloc[:,0] = maf_var_spl_1alt.iloc[:,0].str.replace('chr','')
    maf_var_spl_1alt["VEP"] = ""
    maf_var_spl_1alt = maf_var_spl_1alt[['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER',
                                         'VEP', 'AD_alt', 'DP', 'VAF_alt', 'VAF_ref', 'ALT_num',
                                         'var_type','GT', 'AD', 'INFO', 'FORMAT']]

    ### 7. Save file
    # Get file name without extension form initial path and create temp output folder+file
    file = input_maf.split("/")[-1].split(".")[0]
    temp_file = output_folder + file + ".maf"
    # Save
    maf_var_spl_1alt.to_csv(temp_file, sep="\t", index=False, compression='gzip')

if __name__ == "__main__":
    filter_chgenes()
