import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib import collections as mc
from matplotlib import gridspec
from collections import defaultdict
from scipy.stats import norm
import matplotlib as mpl
import os

import conf as conf

np.random.seed(42)


def get_position(row):
    
    try:
        v = int("".join(row["aachange"][1:-1]))
        return v
    except:
        return -1
    
    
def get_PFAMs_per_transcript(PFAM_files, PFAM_info, transcript):
    df_pfam = pd.read_csv(PFAM_files, sep="\t", names=["ENSEMBL_GENE", "ENSEMBL_TRANSCRIPT", "START", "END", "DOMAIN"])
    df_names = pd.read_csv(PFAM_info, sep="\t", names=["DOMAIN", "CLAN", "CLAN_NAME", "DOMAIN_NAME", "Long Name"])

    # Get domains
    df_pfam_gene = df_pfam[(df_pfam["ENSEMBL_TRANSCRIPT"] == transcript)]
    df_pfam_gene = df_pfam_gene[["ENSEMBL_TRANSCRIPT", "START", "END", "DOMAIN"]].drop_duplicates()
    df_pfam_gene = pd.merge(df_pfam_gene, df_names[["DOMAIN", "DOMAIN_NAME"]].drop_duplicates(), how="left")
    df_pfam_gene["POS"] = df_pfam_gene.apply(lambda row: row["START"] + ((row["END"] - row["START"]) // 2), axis=1)
    df_pfam_gene["SIZE"] = df_pfam_gene.apply(lambda row: row["END"] - row["START"] + 1, axis=1)
    df_pfam_gene["Color"] = "#998ec3"

    return df_pfam_gene



def get_positions_in_CDS(transcript, path_coord):
    df = pd.read_csv(path_coord, sep='\t', low_memory=False,
                     names=['gene', 'gene_symbol', 'prot', 'chr', 's', 'e', 'aa', 'cds', 'genpos',
                            'strand', 'transcript'])

    toappend = []
    strand = ''
    for i, row in df[df['transcript'] == transcript].sort_values(by='s').iterrows():
        toappend.extend([i for i in range(row['s'], row['e'] + 1)])
        strand = row['strand']
    if strand == -1:
        toappend = toappend[::-1]

    return toappend


def load_saturation_CH(path, gene, shap_corrected):

    df = path.copy()
    df.drop_duplicates(inplace=True)
    df["Protein_position"] = df.apply(lambda row: get_position(row), axis=1)

    # aggregate PTMs
    df["PTM"] = df.apply(lambda r: max(r['Phosphorylation'],
                                           r['Acetylation'],
                                           r['Methylation'],
                                           r['Ubiquitination'],
                                           r['Regulatory_Site']), axis=1)
    # SHAP-corrected feature values
    if shap_corrected==True:
        for c in df.columns:
            if c.startswith('shap_'):
                df[c[5:]] = df[c].values * df[c[5:]].values  # multiply each feature by its corresponding SHAP
    # summarize to codon level
    df = df[~df['aachange'].isnull()]
    if len(df) > 0:
        df["AA"] = df.apply(lambda row: row["aachange"][0], axis=1)
        df3 = df.groupby(["AA", "Protein_position", "ENSEMBL_TRANSCRIPT"],
                         as_index=False).agg(
                            {
                                "boostDM_score": np.nanmax,
                                "boostDM_class": np.any,
                                "HotMaps": np.nanmax,
                                "smRegions": np.nanmax,
                                "CLUSTL": np.nanmax,
                                "csqn_type_missense": np.nanmean,
                                "csqn_type_nonsense": np.nanmean,
                                "csqn_type_splicing": np.nanmean,
                                "PhyloP": np.nanmean,
                                "PTM": np.nanmax,
                                "nmd": np.nanmax
                            })
        df3["gene"] = gene
               
        
        df3.sort_values(by='Protein_position', ascending=True, inplace=True)
        df3.reset_index(inplace=True)
        return df3
    print(f"file {path_file} not found...")
    return pd.DataFrame([])

def load_saturation_cancer(path, gene, shap_corrected):

    df = path.copy()
    df.drop_duplicates(inplace=True)
    df["Protein_position"] = df.apply(lambda row: get_position(row), axis=1)

    # aggregate PTMs
    df["PTM"] = df.apply(lambda r: max(r['Phosphorylation'],
                                           r['Acetylation'],
                                           r['Methylation'],
                                           r['Ubiquitination'],
                                           r['Regulatory_Site']), axis=1)

    # summarize to codon level
    df = df[~df['aachange'].isnull()]
    if len(df) > 0:
        df["AA"] = df.apply(lambda row: row["aachange"][0], axis=1)
        df3 = df.groupby(["AA", "Protein_position", "ENSEMBL_TRANSCRIPT"],
                         as_index=False).agg(
                            {
                                "boostDM_score": np.nanmax,
                                "boostDM_class": np.any,
                                "HotMaps": np.nanmax,
                                "smRegions": np.nanmax,
                                "CLUSTL": np.nanmax,
                                "csqn_type_missense": np.nanmean,
                                "csqn_type_nonsense": np.nanmean,
                                "csqn_type_splicing": np.nanmean,
                                "PhyloP": np.nanmean,
                                "PTM": np.nanmax,
                                "nmd": np.nanmax
                            })
        df3["gene"] = gene
        
        df3.sort_values(by='Protein_position', ascending=True, inplace=True)
        df3.reset_index(inplace=True)
        return df3
    print(f"file {path_file} not found...")
    return pd.DataFrame([])


def tracked_blueprint(gene, model,  CH_mat, cancer_mat, df_CH, df_cancer, cancer_name, show=True):

    wanted_df = df_cancer[(df_cancer['gene'] == gene)]

    for transcript, gene in wanted_df[["ENSEMBL_TRANSCRIPT",
                                              "gene"]].drop_duplicates().values:

        # get PFAM domains and subset the mutation data
        subset_data_pfam = get_PFAMs_per_transcript(conf.PFAM_files, conf.PFAM_info, transcript)
        subset_CH_muts = df_CH[
            (df_CH["ENSEMBL_TRANSCRIPT"] == transcript)].sort_values(by='Protein_position',
                                                                         ascending=True)
        subset_cancer_muts = df_cancer[
            (df_cancer["ENSEMBL_TRANSCRIPT"] == transcript)].sort_values(by='Protein_position',
                                                                         ascending=True)

        # define figure layout
        fig = plt.figure(figsize=(8, 3))

        # grid layout
        gs = gridspec.GridSpec(15, 13, figure=fig)

        border = -(len(conf.features_names) + 1)  # bottom limit for blueprint scatter
        
        # MAin plot
        ax1 = plt.subplot(gs[:3, :12])
        #Blueprint
        ax3 = plt.subplot(gs[3, :12], sharex=ax1)

        axes = []
        for i, track in enumerate(conf.features_names):
            axes.append(plt.subplot(gs[border+i+1, :12], sharex=ax1))

        plot_codon_bands(subset_data_pfam, CH_mat, cancer_mat, cancer_name, ax1, ax3)

        for i, track in enumerate(conf.features_names):

            axes[i].plot(subset_CH_muts[track].values, color='tab:red', lw=0.5, alpha=0.5)
            axes[i].plot(subset_cancer_muts[track].values, color='#0000ff', lw=0.5, alpha=0.5)
            axes[i].spines['bottom'].set_visible(False)
            axes[i].spines['left'].set_linewidth(1)
            axes[i].spines['right'].set_visible(False)
            axes[i].spines['top'].set_visible(False)
            axes[i].set_yticks([])
            axes[i].get_xaxis().set_visible(False)
            axes[i].set_ylabel(conf.features_names[track], rotation=0, labelpad=7, fontsize=6,
                               horizontalalignment='right', verticalalignment='center')

        ax1.set_title(f'{gene}', fontsize=10)
        if show:
            plt.show()
        plt.close(fig)
        

def plot_codon_bands(df_pfam_gene, CH_mat, cancer_mat, cancer_name, ax_0, ax_2):

#     ax_0.set_xticks(np.linspace(0, 1, 3))
    ax_0.spines['bottom'].set_visible(False)
    ax_0.spines['left'].set_linewidth(1)
    ax_0.spines['right'].set_visible(False)
    ax_0.spines['top'].set_visible(False)
    CH_mat["Protein_position"] = CH_mat.apply(lambda row: get_position(row), axis=1)
    cancer_mat["Protein_position"] = cancer_mat.apply(lambda row: get_position(row), axis=1)
   
    # set equivalent coordinates for the three possible mutations
    ax_0.plot(CH_mat[CH_mat['boostDM_class']==True]['Protein_position'], np.full_like(CH_mat[CH_mat['boostDM_class']==True]['Protein_position'],2), '|', markeredgewidth=0.2, markersize=6, color='tab:red')
    ax_0.plot(cancer_mat[cancer_mat['boostDM_class']==True]['Protein_position'], np.full_like(cancer_mat[cancer_mat['boostDM_class']==True]['Protein_position'],1), '|k', markeredgewidth=0.2,  markersize=6, color='#0000ff')
    ax_0.set_yticks(np.arange(3)[1:])
    ax_0.set_yticklabels(['Myeloid Cancer', 'CH'], fontsize = 6)
    ax_0.set_xticks(np.append(np.arange(0, max(CH_mat['Protein_position']), 100)[:-1], max(CH_mat['Protein_position'])))
    ax_0.set_xticklabels(np.append(np.arange(0, max(CH_mat['Protein_position']), 100)[:-1], max(CH_mat['Protein_position'])), fontsize = 0)
    ax_0.set_ylim(.5, 2.5)
    ax_0.set_xlim(0, max(CH_mat['Protein_position']))
    ax_2.set_ylim(0, 1)
    
    for i, r in df_pfam_gene.iterrows():
        start_base = r['START']
        size_base = r['SIZE']
        rect = patches.Rectangle(xy=(start_base, 0), width=size_base, height=5, color=r["Color"], alpha=0.5, zorder=2)
        ax_2.annotate(text=r["DOMAIN_NAME"], xy=(start_base + 1, 0.3), fontsize=5)
        ax_2.add_patch(rect)

    ax_2.set_xticklabels(np.append(np.arange(0, max(CH_mat['Protein_position']), 100)[:-1], max(CH_mat['Protein_position'])), fontsize = 6)
    ax_2.set_yticks([])

    
def plot_gene_full_nucleotide_invert(subset_data_pfam, cancer_type, df, transcript, path_coord,\
                              ax0, all_possible=False):

    # remove those mutations not falling in CDS:
    df = df[df['AA'] != 'n']

    # Configure the axis
    ax0.set_title('Observed mutations in Myeloid Cancer (n='+str(df["number_observed_muts"].sum())+')', y=-0.25)
    ax0.set_ylabel("mutation count")

    ax0.spines['bottom'].set_visible(False)
    ax0.spines['left'].set_linewidth(1)
    ax0.spines['right'].set_visible(False)
    ax0.spines['top'].set_linewidth(1)
    ax0.tick_params(axis='y', labelsize=6, pad=0.25, width=0.25, length=1.5)
    ax0.set_xticks([])

    # set equivalent coordinates for the three possible mutations
    set_coordinates = get_positions_in_CDS(transcript, path_coord)


    # we need to get the set of equivalent coordinates per gene
    equivalent_coordinates = {coord: i for i, coord in enumerate(set_coordinates)}
    vals_coord = list(equivalent_coordinates.values())
    axs = [ax0]
    for ax in axs:
        ax.set_xlim(np.min(vals_coord), np.max(vals_coord))

    # plot observed mutations
    pos_list = df["pos"].tolist()
    ys = df["number_observed_muts"].values

    d = df["boostDM_score"].values

    coordinates_mutations = []

    passenger_x = []
    passenger_y = []
    passenger_color = []

    driver_x = []
    driver_y = []
    driver_color = []

    # for each of the positions
    for i, p in enumerate(pos_list):
        if ys[i] > 0:

            coordinates_mutations.append([(equivalent_coordinates[p], 0),
                                          (equivalent_coordinates[p], ys[i] - 0.1)])

            if d[i] < 0.5:

                passenger_x.append(equivalent_coordinates[p])
                if all_possible:
                    passenger_y.append(d[i])
                else:
                    passenger_y.append(ys[i])
                passenger_color.append('#636363')

            else:
                driver_x.append(equivalent_coordinates[p])
                if all_possible:
                    driver_y.append(d[i])
                else:
                    driver_y.append(ys[i])
                driver_color.append('#ac0f0f')

    lc = mc.LineCollection(coordinates_mutations, colors='black', linewidths=1, alpha=0.3)
    ax0.add_collection(lc)

    size = 12
    ax0.scatter(passenger_x, passenger_y, s=size, c=passenger_color, alpha=0.7, label='non-driver')
    ax0.scatter(driver_x, driver_y, s=size, c=driver_color, alpha=0.7, label='driver')
    leg = ax0.legend()
    leg.get_frame().set_linewidth(0.0)
    
    ax0.invert_yaxis()
    ax0.xaxis.tick_top()



def plot_gene_full_nucleotide_simple(subset_data_pfam, df, transcript, path_coord,
                              ax0, all_possible=False):

    # remove those mutations not falling in CDS:
    df = df[df['AA'] != 'n']

    # Configure the axis
    ax0.set_title('Observed mutations in CH (n='+str(df["number_observed_muts"].sum())+')')
    ax0.set_ylabel("mutation count")

    ax0.spines['bottom'].set_visible(1)
    ax0.spines['left'].set_linewidth(1)
    ax0.spines['right'].set_visible(False)
    ax0.spines['top'].set_linewidth(False)
    ax0.tick_params(axis='y', labelsize=6, pad=0.25, width=0.25, length=1.5)

    # set equivalent coordinates for the three possible mutations
    set_coordinates = get_positions_in_CDS(transcript, path_coord)


    # we need to get the set of equivalent coordinates per gene
    equivalent_coordinates = {coord: i for i, coord in enumerate(set_coordinates)}
    vals_coord = list(equivalent_coordinates.values())
    axs = [ax0]
    for ax in axs:
        ax.set_xlim(np.min(vals_coord), np.max(vals_coord))

    # plot observed mutations
    pos_list = df["pos"].tolist()
    ys = df["number_observed_muts"].values

    d = df["boostDM_score"].values

    coordinates_mutations = []

    passenger_x = []
    passenger_y = []
    passenger_color = []

    driver_x = []
    driver_y = []
    driver_color = []

    # for each of the positions
    for i, p in enumerate(pos_list):
        if ys[i] > 0:

            coordinates_mutations.append([(equivalent_coordinates[p], 0),
                                          (equivalent_coordinates[p], ys[i] - 0.1)])

            if d[i] < 0.5:

                passenger_x.append(equivalent_coordinates[p])
                if all_possible:
                    passenger_y.append(d[i])
                else:
                    passenger_y.append(ys[i])
                passenger_color.append('#636363')

            else:
                driver_x.append(equivalent_coordinates[p])
                if all_possible:
                    driver_y.append(d[i])
                else:
                    driver_y.append(ys[i])
                driver_color.append('#ac0f0f')

    lc = mc.LineCollection(coordinates_mutations, colors='black', linewidths=1, alpha=0.3)
    ax0.add_collection(lc)

    size = 12
    ax0.scatter(passenger_x, passenger_y, s=size, c=passenger_color, alpha=0.7, label='non-driver')
    ax0.scatter(driver_x, driver_y, s=size, c=driver_color, alpha=0.7, label='driver')
    leg = ax0.legend(loc='upper left')
    leg.get_frame().set_linewidth(0.0)
    
    
def plot_observed_distribution_figure_3_simple(gene, CH_mat, cancer_mat, cancer_type,
                                        data_ch, count_driver_ch, count_driver_unique_ch, count_total_ch,
                                        count_total_unique_ch, score_values_ch, count_missense_ch, count_nonsense_ch,
                                        count_synonymous_ch,
                                        data_cancer, count_driver_cancer, count_driver_unique_cancer,
                                        count_total_cancer, count_total_unique_cancer, score_values_cancer,
                                        count_missense_cancer, count_nonsense_cancer, count_synonymous_cancer,
                                    dpi=200):

    conf.config_params(font_size=6)

    mpl.rcParams['figure.dpi']= dpi
    wanted_df = data_ch[data_ch['gene'] == gene]

    for transcript, gene in wanted_df[["ENSEMBL_TRANSCRIPT","gene"]].drop_duplicates().values:

        # get PFAM domains and subset the mutation data
        subset_data_pfam = get_PFAMs_per_transcript(conf.PFAM_files, conf.PFAM_info, transcript)
        subset_data_muts_ch = data_ch[data_ch["ENSEMBL_TRANSCRIPT"] == transcript]
        subset_data_muts_cancer = data_cancer[data_cancer["ENSEMBL_TRANSCRIPT"] == transcript]


        # define figure layout
        fig = plt.figure(figsize=(5, 4))
        fig.suptitle(gene, y=0.95)
        gs = gridspec.GridSpec(19, 3, figure=fig)
                
        ax1 = plt.subplot(gs[1:6,:])
        ax4 = plt.subplot(gs[6:9, :])
        ax5 = plt.subplot(gs[9,:])
        ax7 = plt.subplot(gs[11:16, :],sharex=ax1)

        plot_gene_full_nucleotide_simple(subset_data_pfam, subset_data_muts_ch, transcript, conf.path_coord, ax1)
        plot_codon_bands_joanen_simple(subset_data_pfam, CH_mat, cancer_type, cancer_mat, ax4, ax5)
        plot_gene_full_nucleotide_invert(subset_data_pfam, cancer_type, subset_data_muts_cancer, transcript, conf.path_coord, ax7)

        plt.show()

        
def plot_codon_bands_joanen_simple(df_pfam_gene, CH_mat, cancer_type, cancer_mat, ax_0, ax_2):

    ax_0.spines['bottom'].set_visible(False)
    ax_0.spines['left'].set_linewidth(1)
    ax_0.spines['right'].set_visible(False)
    ax_0.spines['top'].set_visible(False)
    CH_mat["Protein_position"] = CH_mat.apply(lambda row: get_position(row), axis=1)
    cancer_mat["Protein_position"] = cancer_mat.apply(lambda row: get_position(row), axis=1)
   
    # set equivalent coordinates for the three possible mutations
    ax_0.plot(CH_mat[CH_mat['boostDM_class']==True]['Protein_position'], np.full_like(CH_mat[CH_mat['boostDM_class']==True]['Protein_position'],2), '|', markeredgewidth=0.1, markersize=6, color='#ac0f0f')
    ax_0.plot(cancer_mat[cancer_mat['boostDM_class']==True]['Protein_position'], np.full_like(cancer_mat[cancer_mat['boostDM_class']==True]['Protein_position'],1), '|k', markeredgewidth=0.1,  markersize=6, color='#0000ff')
    ax_0.set_yticks(np.arange(3)[1:])
    ax_0.set_yticklabels(['Myeloid Cancer', 'CH'], fontsize = 6)
    ax_0.set_xticks([])

    ax_0.set_ylim(.5, 2.5)
    ax_0.set_xlim(2, max(CH_mat['Protein_position']))
    
    ax_2.set_ylim(0, 1)

    for i, r in df_pfam_gene.iterrows():
        start_base = r['START']
        size_base = r['SIZE']
        rect = patches.Rectangle(xy=(start_base, 0), width=size_base, height=5, color=r["Color"], alpha=0.5, zorder=2)
        ax_2.annotate(text=r["DOMAIN_NAME"], xy=(start_base + 1, 0.3), fontsize=5)
        ax_2.add_patch(rect)

    ax_2.set_xticks(np.append(np.arange(0, max(CH_mat['Protein_position']), 100)[:-1], max(CH_mat['Protein_position'])))
    ax_2.set_xticklabels(np.append(np.arange(0, max(CH_mat['Protein_position']), 100)[:-1], max(CH_mat['Protein_position'])), fontsize = 6)
    ax_2.set_xlim(2, max(CH_mat['Protein_position']))
    ax_2.set_yticks([])