import os
import matplotlib as mpl
import matplotlib.pyplot as plt


# check these paths

output_boostdm = '../../Paper_data/BoostDM-CH/'
output_boostdm_cancer = '../../Paper_data/BoostDM-cancer/'
intogen_cancer_path = '../../Paper_data/Intogen-Cancer/'

#PFAM info
PFAM_files = '../../Paper_data/pfam/pfam_biomart.tsv.gz'
PFAM_info = '../../Paper_data/pfam/pfam_names.info.csv'
# CDS coordinates
path_coord = '../../Paper_data/pfam/cds_biomart.tsv'

# #new
CANCER_COHORTS_PATH = os.path.join(intogen_cancer_path, 'cohorts.tsv')
CANCER_DRIVERS_PATH = os.path.join(intogen_cancer_path, 'drivers.tsv')
ONCOTREE_PATH = os.path.join(intogen_cancer_path, 'oncotree-tree.tsv')


features = ["shap_CLUSTL_SCORE","shap_CLUSTL_cat_1","shap_HotMaps_cat_1",
            "shap_smRegions_cat_1","shap_csqn_type_missense","shap_csqn_type_nonsense",
            "shap_csqn_type_splicing","shap_PhyloP", "shap_Acetylation","shap_Phosphorylation",
            "shap_Regulatory_Site","shap_Ubiquitination","shap_Methylation","shap_nmd"]

# naming conventions
features_names = {'HotMaps': '3D cluster',
         'CLUSTL': 'Linear cluster',
         'smRegions': 'PFAM Domain',
         'PhyloP': 'Conservation',
         'PTM': 'PTM',
         'nmd': 'NMD',
         'csqn_type_missense': 'Missense',
         'csqn_type_nonsense': 'Nonsense',
         'csqn_type_splicing': 'Splicing'}


names_features = {
    'shap_nmd': 'NMD',
    'shap_CLUSTL_int': 'Linear cluster',
    'shap_HotMaps_cat_1': '3D cluster',
    'shap_smRegions_cat_1': 'PFAM domain',
    'shap_PTM_int': 'PTM',
    'shap_PhyloP': 'Conservation',
    'shap_csqn_type_missense': 'Missense',
    'shap_csqn_type_splicing': 'Splicing',
    'shap_csqn_type_nonsense': 'Nonsense'}


def config_params(font_size=12):

    mpl.rcParams.update(mpl.rcParamsDefault)
    plt.rcParams['font.sans-serif'] = ['arial']
    plt.rcParams['font.size'] = font_size
    plt.rcParams['font.family'] = ['sans-serif']
    plt.rcParams['svg.fonttype'] = 'none'
    plt.rcParams['mathtext.fontset'] = 'custom'
    plt.rcParams['mathtext.cal'] = 'arial'
    plt.rcParams['mathtext.rm'] = 'arial'
