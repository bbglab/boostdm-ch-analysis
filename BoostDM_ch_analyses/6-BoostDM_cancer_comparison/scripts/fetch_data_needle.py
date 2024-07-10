import sys
import os

import seaborn as sns
sns.set_style("whitegrid", {'axes.grid' : False})
import pandas as pd
import numpy as np

import conf
import conf
import oncotree_notebook
tree = oncotree_notebook.Oncotree()


# Take observed mutations and merges the prediciton from quality models
def create_observed_dataset(prediction_path, gene, cohort, obs_mut):
    
    sat_pred = pd.read_csv(prediction_path, sep='\t')
    sat_pred['chr'] = sat_pred['chr'].astype(str)
    sat_pred['pos'] = sat_pred['pos'].astype(int)
    obs_mut = obs_mut[obs_mut['gene'] == gene]
    df = obs_mut.merge(sat_pred, on=['gene', 'chr', 'pos', 'alt', 'aachange'], how='left')
    df.rename(columns={"COHORT": "cancer_type"}, inplace=True)
    df = df[~df['boostDM_class'].isnull()]
    return df

# Take observed mutations from cancer
def create_observed_dataset_cancer(prediction_path, gene, cohort, obs_mut):

    sat_pred = pd.read_csv(prediction_path, sep='\t')
    sat_pred['chr'] = sat_pred['chr'].astype(str)
    sat_pred['pos'] = sat_pred['pos'].astype(int)
    # ! SDM CHANGE
    cohort_list = tree.get_cohorts(cohort)
    cohort_list = '|'.join(x for x in cohort_list)
    obs_mut = obs_mut[(obs_mut['gene'] == gene) & (obs_mut['COHORT'].str.contains(cohort_list))]
#     obs_mut = obs_mut[obs_mut['gene'] == gene]
    df = obs_mut.merge(sat_pred, on=['gene', 'chr', 'pos', 'alt', 'aachange'], how='left')
    df.rename(columns={"COHORT": "cancer_type"}, inplace=True)
    df = df[~df['boostDM_class'].isnull()]
    return df


def get_position(row):
    
    try:
        v = int("".join(row["aachange"][1:-1]))
        return v
    except:
        return -1


def get_plot_data_joanen(data):

    data["Protein_position"] = data.apply(lambda row: get_position(row), axis=1)
    data["AA"] = data.apply(lambda row: row["aachange"][0], axis=1)
    data['ID'] = data.apply(lambda x: '{}_{}'.format(x['pos'], x['alt']), axis=1)

    score_values = data['boostDM_score'].tolist()
    count_driver_unique = len(set(data[data['boostDM_class']]['ID']))  # predicted to be drivers (unique)
    count_driver = len(data[data['boostDM_class']]['ID'])  # predicted to be drivers
    count_total = len(data['ID'])  # total number of mutations
    count_total_unique = len(set(data['ID']))  # total mutations (unique)
    count_missense = len(data[data['impact']=='Missense']['ID'])
    count_nonsense = len(data[data['impact']=='Nonsense']['ID'])
    count_synonymous = len(data[data['impact']=='Synonymous']['ID'])

    data = data.groupby(["ID", "pos", "AA", "Protein_position", "gene", "ENSEMBL_TRANSCRIPT", "boostDM_score", "boostDM_class", "impact"], as_index=False).agg({"sampleID": "count"})
    data.rename(columns={"sampleID": "number_observed_muts"}, inplace=True)

    return data, count_driver, count_driver_unique, count_total, count_total_unique, score_values,\
count_missense, count_nonsense, count_synonymous 
