import os
import pandas as pd
import numpy as np

from gene_utils import *
from figure_utils import *
from constants import *

import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

mpl.style.use(os.path.join(ASSETS_DIR, 'stylesheet.mplstyle'))
    
def make_fig2_class_enrichment_legend():
    """
    Generate a figure legend for barplot colors in Figure 2b
    """
    handles = ['In class', 'Out of class']
    plt.figure(figsize=(15 * mm, 7 * mm))
    plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
    handles = [
        Patch(facecolor='tab:orange', edgecolor='white', label='In class'),
        Patch(facecolor='tab:blue', edgecolor='white', label='Out of class')
    ]
    leg = plt.legend(handles = handles, ncol=1, 
                     handlelength=1, handleheight=1,
                     handletextpad=0.5, columnspacing=1, loc='lower left', prop={'size': ANNOT_SIZE})
    plt.gca().spines[['left', 'right', 'top', 'bottom']].set_visible(False)
    plt.gca().set_xticks([])
    plt.gca().set_yticks([])
    plt.savefig(os.path.join(FIGURE_DIR, 'fig2_gene_class_legend.pdf'))
    
def make_fig2_vertical_lineage_legend():
    """
    Generate a figure legend for lineages of scatterplot points across Figure 2
    """
    lineages_to_show = ['Ampulla of Vater', 'Breast', 'Biliary Tract', 'Colorectal', 'CNS/Brain', 
                        'Esophagus/Stomach', 'Ovary/Fallopian Tube', 'Pancreas', 'Prostate', 'Uterus']

    plt.figure(figsize=(25 * mm, 35 * mm))
    plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
    handles = []
    # handles.append(Patch(facecolor='white', edgecolor='white', label='Geneset'))
    for gs in lineages_to_show:
        c = 'black'
        handle = Line2D([], [], color=lineage_cmap[gs], marker = 'o', linestyle='None', markersize=4, markeredgewidth=0, label = gs)
        handles.append(handle)    
    handles.append(Line2D([], [], color='black', marker = 's', linestyle='None', markersize=4, label = 'NextGen', markeredgewidth=0))
    handles.append(Line2D([], [], color='black', marker = 'o', linestyle='None', markersize=4, label = 'Adherent', markeredgewidth=0))
    leg = plt.legend(handles = handles, ncol=1, handletextpad=0, columnspacing=1, loc='lower left', prop={'size': ANNOT_SIZE})
    plt.gca().spines[['left', 'right', 'top', 'bottom']].set_visible(False)
    plt.gca().set_xticks([])
    plt.gca().set_yticks([])
    plt.savefig(os.path.join(FIGURE_DIR, 'fig2_lineage_legend.pdf'))
    
def make_fig4_wide_lineage_legend():
    """
    Generate a figure legend for lineages of scatterplot points across Figure 4
    """
    lineages_to_show = ['Ampulla of Vater', 'Breast', 'Biliary Tract', 'Colorectal', 'Esophagus/Stomach', 'Ovary/Fallopian Tube', 'Pancreas', 'Prostate', 'Uterus']

    plt.figure(figsize=(55 * mm, 10 * mm))
    plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
    handles = []
    # handles.append(Patch(facecolor='white', edgecolor='white', label='Geneset'))
    for gs in lineages_to_show:
        c = 'black'
        handle = Line2D([], [], color=lineage_cmap[gs], marker = 'o', linestyle='None', markersize=4, markeredgewidth=0, label = gs)
        handles.append(handle)    
    leg = plt.legend(handles = handles, ncol=3, handletextpad=0, columnspacing=1, loc='lower left', prop={'size': ANNOT_SIZE})
    plt.gca().spines[['left', 'right', 'top', 'bottom']].set_visible(False)
    plt.gca().set_xticks([])
    plt.gca().set_yticks([])
    plt.savefig(os.path.join(FIGURE_DIR, 'fig4_lineage_legend.pdf'))
    
def make_fig5_geneset_legend():
    """
    Generate a figure legend for gene sets corresponding to scatterplot points across Figure 5
    """
    genesets_to_show = ['Integrin', 'Actin Regulation', 'Adherens Junction', 'Tight Junction', 'Lipid Metabolism', 'Wnt Signaling']

    plt.figure(figsize=(80 * mm, 7 * mm))
    plt.subplots_adjust(left=0.05, right=0.85, top=0.6, bottom=0)
    handles = []
    for gs in genesets_to_show:
        c = 'black'
        handle = Line2D([], [], color=mp_gs_pal[gs], marker = 'o', linestyle='None', markersize=4, markeredgewidth=0, label = gs)
        handles.append(handle)    
    handles.append(Line2D([], [], color='tab:gray', marker = 'o', linestyle='None', markersize=4, label = 'Other', markeredgewidth=0, alpha=0.8))
    handles.append(Line2D([], [], color='tab:gray', marker = 'o', linestyle='None', markersize=4, label = 'Nondependency', markeredgewidth=0, alpha=0.3))
    leg = plt.legend(handles = handles, ncol=4, handletextpad=0, columnspacing=1, loc='lower left', prop={'size': ANNOT_SIZE})
    plt.gca().spines[['left', 'right', 'top', 'bottom']].set_visible(False)
    plt.gca().set_xticks([])
    plt.gca().set_yticks([])
    plt.savefig(os.path.join(FIGURE_DIR, 'fig5_minipool_geneset_legend.pdf'))

def main():
    """
    Generate all missing figure legends
    """
    make_fig2_class_enrichment_legend()
    make_fig2_vertical_lineage_legend()
    make_fig4_wide_lineage_legend()
    make_fig5_geneset_legend()
    

if __name__ == "__main__":
    main()