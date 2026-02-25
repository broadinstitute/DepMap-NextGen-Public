import os

import matplotlib.pyplot as plt
import seaborn as sns

INPUT_DIR = os.path.join('data')
PROCESSED_DIR = os.path.join('processed')
FIGURE_DIR = os.path.join('figures')
ASSETS_DIR = os.path.join('src', 'assets')

# conversions from inches
cm = 1/2.54
mm = cm/10
pt = 1/72

# font sizes (in cases where we manually override the stylesheet)
TITLE_SIZE = 7
LABEL_SIZE = 6
TICK_SIZE = 5
ANNOT_SIZE = 5
POINT_SIZE = 8

# annotation arrow styles
annot_arrow_props = dict(arrowstyle='-', color='black', alpha=1, linewidth=0.25, shrinkA=1, shrinkB=1)

# lineage renaming
lineage_replacement = {'Bowel': 'Colorectal'}

# colormaps
lineage_cmap = {
    'Esophagus/Stomach': '#018e00',
    'Pancreas': '#521b92',
    'Biliary Tract': '#FAB19B',
    'CNS/Brain': '#0095ff',
    'Breast': '#005492',
    'Bowel': '#ff9200',
    'Colorectal': '#ff9200',
    'Prostate': '#931100',
    'Ovary/Fallopian Tube': '#ff2600',
    'Ampulla of Vater': '#011892',
    'Uterus': '#918f00',
    'Lung': '#7f7f7f', 
    'Head and Neck': '#7f7f7f',
    'Other': '#7f7f7f'
}

oncoprint_palette = {
    'Amplification': plt.cm.coolwarm(256), #'tab:red', 
    'Deletion': plt.cm.coolwarm(0), #'tab:blue', 
    'Neutral Copy': 'white'
}

dependency_class_palette = {
    'Growth Suppressor': '#812afa',
    'Non-dependency': '#00204DFF',
    'Weakly Selective': 'tab:gray',
    'Strongly Selective': '#145A32',
    'High Variance': '#FF6F00B2',
    'Pan-dependency': '#C71000B2'
}

addiction_class_cmap = {
    'Oncogenic TF': '#004D40', 
    'Oncogene': '#1E88E5', 
    'Transcription factor': '#FFC107', 
    'TF': '#FFC107', 
    'Other': 'tab:gray'
}

onc_gof_tsg_lof_cmap = {
    'Insignificant': '#cdd7d6', 
    'Driven by alteration': '#048ba8',
    'Driven by lineage': '#898c8c'
}

depmap_set_palette = {
    'Cas9 Adherent': sns.color_palette("Dark2")[0], 
    'Humagne Adherent': sns.color_palette("Dark2")[5], 
    'Cas9 2D': sns.color_palette("Dark2")[0], 
    'Humagne 2D': sns.color_palette("Dark2")[5], 
    'Cas9 Traditional': sns.color_palette("Dark2")[0], 
    'Humagne Traditional': sns.color_palette("Dark2")[5], 
    'Organoids': sns.color_palette("Dark2")[2], 
    'Neurospheres': sns.color_palette("tab10")[9],
    'Spheroids': sns.color_palette("tab10")[9]
}

sample_set_palette = {
    'TCGA+': 'tab:orange', 
    'DepMap Organoid': 'tab:blue', 
    'DepMap Adherent': 'tab:green',
    'DepMap 2D': 'tab:green'
}

cluster_cmap = {
    '1': '#06d6a0', # green 
    'SLE': '#06d6a0',
    '2': '#ef476f', # red 
    'Lipid Metabolism': '#ef476f',
    '3': '#f78c6b', # orange 
    'OxPhos': '#f78c6b',
    '4': '#ffd166', # yellow 
    'Cell Cycle': '#ffd166',
    '5': '#ac44ff', # purple 
    'mTORC1': '#ac44ff',
    '6': '#073b4c', # slate 
    'Adhesion/Cytoskeleton': '#073b4c',
    '7': '#118ab2', # blue 
    'UV': '#118ab2',
    'Other': 'tab:gray',
}

mp_gs_pal = {
    'Integrin': '#0cb2af',
    'Actin Regulation': '#a1c65d',
    'Actin Nucleation': 'tab:gray',
    'Adherens Junction': '#f29222',
    'Tight Junction': '#e95e50',
    'Lipid Metabolism': '#f18aad',
    'Wnt Signaling': '#936fac',
    'Other': 'tab:gray'
}

viability_palette = {
    'Plastic': '#f58f29',
    'Dome-R10': '#a4b0f5',
    'Dome-OPAC': '#4464ad'
}

highlight_color = '#f87060'