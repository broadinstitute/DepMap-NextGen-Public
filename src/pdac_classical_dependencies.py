import os
import pandas as pd
import numpy as np
from figure_utils import *
from constants import *
from data_utils import *

import seaborn as sns
import matplotlib.pyplot as plt

mpl.style.use(os.path.join(ASSETS_DIR, 'stylesheet.mplstyle'))

MANUAL_ANNOTATE = True

# figure 4d

def plot_pdac_volcano(fdr_threshold=0.05, wnt_color='darkturquoise', highlight_color='#f87060', manual_annotate=MANUAL_ANNOTATE):
    """
    Plot correlations for dependencies against the expression of the PDAC-Classical program (from Gavish et al., 2023), and
    show the running enrichment score for Wnt related genes in parallel

    Args:
        fdr_threshold (float, optional): Minimum significance value. Defaults to 0.05.
        wnt_color (str, optional): String color to use for highlighting a geneset. Defaults to 'darkturquoise'.
        highlight_color (str, optional): String color to use for manual highlighting of points. Defaults to '#f87060'.
        manual_annotate (bool, optional): Flag for whether to manually annotate genes. Defaults to MANUAL_ANNOTATE.
    """
    pdac_classical_dependency_df = load_file('pdac_classical_dependency_table.csv', local_dir=PROCESSED_DIR, index_col=0)
    pdac_classical_dependency_df['p_transform'] = -np.log10(pdac_classical_dependency_df['NextGenPValue'])
    pdac_classical_dependency_df['label'] = pdac_classical_dependency_df.index.str.split(' ').str[0]
    
    fig, axs = plt.subplots(3, 1, figsize=(50 * mm, 60 * mm), height_ratios=[20, 4, 2], gridspec_kw={'hspace':0})
    plt.subplots_adjust(right=0.95, top=0.88, bottom=0.15)

    xlabel='NextGenPearsonR'
    ylabel='p_transform'
    
    # plot the volcano
    sns.scatterplot(
        pdac_classical_dependency_df[~pdac_classical_dependency_df['Wnt Signaling']].sort_values(['p_transform']),
        x=xlabel, y=ylabel,
        color='lightgray',
        s=4,
        alpha=0.6, edgecolor='black', linewidth=0.25,
        ax=axs[0], legend=False
    )
    sns.scatterplot(
        pdac_classical_dependency_df[pdac_classical_dependency_df['Wnt Signaling']].sort_values(['p_transform']),
        x=xlabel, y=ylabel,
        color=wnt_color,
        s=4,
        alpha=0.6, edgecolor='black', linewidth=0.25,
        ax=axs[0], legend=False
    )
    axs[0].set_xlabel('')
    axs[0].set_xticks([])
    axs[0].set_ylabel('-log10(p-value)', fontdict={'size': LABEL_SIZE})
    axs[0].set_yticks(axs[0].get_yticks()[1:-1], labels=axs[0].get_yticklabels()[1:-1], fontdict={'size': TICK_SIZE})

    # add a significance line
    p_fdr_approx = approximate_fdr_in_p_units(pdac_classical_dependency_df, fdr_threshold=fdr_threshold, pcolumn='NextGenPValue', fdr_column='FDR')
    axs[0].axhline(-np.log10(p_fdr_approx), linestyle='dashed', color='tab:red', alpha=0.5, linewidth=0.5)

    # add annotations
    annots = pdac_classical_dependency_df.query('`Wnt Signaling` & `FDR` < 0.05').sort_values('FDR').set_index('label')
    if manual_annotate:
        manually_annotate(axs[0], 'MESD', annots, 
                          xlabel, ylabel, (0.05, 0), ha='left', va='center',
                          arrowprops=dict(arrowstyle='-', color='black', alpha=1, linewidth=0.25, shrinkA=0, shrinkB=1),
                          color=highlight_color)
        manually_annotate(axs[0], 'WLS', annots, 
                          xlabel, ylabel, (0.05, 0.2), ha='left', va='center',
                          arrowprops=dict(arrowstyle='-', color='black', alpha=1, linewidth=0.25, shrinkA=0, shrinkB=1),
                          color=highlight_color)
        manually_annotate(axs[0], 'TCF7L2', annots, 
                          xlabel, ylabel, (0.05, .1), ha='left', va='center',
                          arrowprops=dict(arrowstyle='-', color='black', alpha=1, linewidth=0.25, shrinkA=0, shrinkB=1),
                          color=highlight_color)
        manually_annotate(axs[0], 'FZD5', annots, 
                          xlabel, ylabel, (0.05, 0), ha='left', va='center',
                          arrowprops=dict(arrowstyle='-', color='black', alpha=1, linewidth=0.25, shrinkA=0, shrinkB=1),
                          color=highlight_color)
        manually_annotate(axs[0], 'LGR4', annots, 
                          xlabel, ylabel, (0.05, 0), ha='left', va='center',
                          arrowprops=dict(arrowstyle='-', color='black', alpha=1, linewidth=0.25, shrinkA=0, shrinkB=1))
        manually_annotate(axs[0], 'SOX9', annots, 
                          xlabel, ylabel, (0.08, .3), ha='left', va='center',
                          arrowprops=dict(arrowstyle='-', color='black', alpha=1, linewidth=0.25, shrinkA=0, shrinkB=1))

    else:
        texts = [axs[0].text(r[xlabel], r[ylabel], i, ha='center', va='center', fontsize=ANNOT_SIZE) for i,r in annots.iterrows()]

    # add a legend
    handles = []
    handles.append(Line2D([], [], color=wnt_color, marker = 'o', linestyle='None', markeredgewidth=0.25, markersize=4, 
                          label='Wnt Signaling', markeredgecolor='black', alpha=0.6))
    axs[0].legend(handles = handles, loc='upper right', prop={'size': ANNOT_SIZE})

    # plot the enrichment line    
    sns.lineplot(
        x=pdac_classical_dependency_df.sort_values(xlabel)[xlabel].values.tolist(),
        y=pdac_classical_dependency_df.sort_values(xlabel)['RunningEnrichmentScore'].values.tolist(),
        color=wnt_color,
        linewidth=0.5,
        ax=axs[1]
    )
    axs[1].axhline(0, linestyle='solid', color='gray', alpha=0.5, linewidth=0.5)
    axs[1].set_ylabel('ES', fontdict={'size': LABEL_SIZE})
    axs[1].set_yticks([0], labels=[0], fontdict={'size': TICK_SIZE})
    axs[1].set_xticks([])
    axs[1].set_xlim(axs[0].get_xlim())

    # plot vertical bars aligned to genes along the enrichment line
    for j, k in pdac_classical_dependency_df[pdac_classical_dependency_df['Wnt Signaling']].iterrows():
        axs[2].axvline(x=k[xlabel], color=wnt_color, linewidth=0.5)
    axs[2].set_yticks([])
    axs[2].set_xlabel('Pearson Correlation', fontdict={'size': LABEL_SIZE})
    axs[2].set_xlim(axs[0].get_xlim())
    axs[2].set_xticks(axs[2].get_xticks()[1:-1], labels=axs[2].get_xticklabels()[1:-1], fontdict={'size': TICK_SIZE})

    axs[0].set_title('PDAC-classical/EED\nAssociated Dependencies', fontdict={'size': TITLE_SIZE})
    plt.savefig(os.path.join(FIGURE_DIR, 'Fig_4d_organoid_pdac_classical_dependency_volcano.pdf'))

# figure 4e

def plot_wnt_regressions(wnt_genes_to_regress):
    """
    Plot scatterplots of example Wnt gene dependencies against PDAC-classical expression

    Args:
        wnt_genes_to_regress (list[str]): A list of genes to make scatteplots for
    """
    example_df = load_file('wnt_example_table.csv', local_dir=PROCESSED_DIR)
    example_df = example_df[example_df['SampleSet'] == "DepMap Organoid"]
    pdac_classical_dependency_df = load_file('pdac_classical_dependency_table.csv', local_dir=PROCESSED_DIR, index_col=0)
    
    # create the subplot series
    fig, axs = plt.subplots(4, 1, figsize=(30 * mm, 60 * mm), sharex=True, sharey=True, gridspec_kw={'hspace': 0.05})
    plt.subplots_adjust(left=0.25, bottom=0.18, top=0.98, right=0.95)
    axs = axs.flatten()

    # makxae pairs of regression plots on top of scatter plots
    for i, g in enumerate(wnt_genes_to_regress):
        example_subset = example_df[example_df['Dependency'] == g]
        
        sns.regplot(
            example_subset,
            x='PDAC-classical',
            y='DependencyGeneEffect',
            scatter=False, color='tab:gray', line_kws={'linewidth': 1},
            ax=axs[i]
        )
        
        sns.scatterplot(
            example_subset,
            x='PDAC-classical',
            y='DependencyGeneEffect',
            hue='OncotreeLineage', palette=lineage_cmap, legend=False, s=4, marker='s',
            ax=axs[i]
        )

        axs[i].text(0.15, 0.06, f'r = {pdac_classical_dependency_df.loc[g, "NextGenPearsonR"]:.2f}', 
                    fontdict={'fontsize': ANNOT_SIZE, 'color': 'black'}, transform=axs[i].transAxes)
        axs[i].set_xlabel('')
        axs[i].set_ylabel('')
        axs[i].set_yticks([-2, -1, 0], [-2, -1, 0], fontdict={'size': TICK_SIZE})
        axs[i].axvline(0, linestyle='dotted', color='black', zorder=-1, linewidth=0.5)
        axs[i].axhline(0, linestyle='dotted', color='black', zorder=-1, linewidth=0.5)

    axs[3].set_xticks([0, 2.5, 5], [0, 2.5, 5], fontdict={'size': TICK_SIZE})
    axs[3].set_xlabel('PDAC-classical/EED Score', x=0.4, fontdict={'size': LABEL_SIZE})

    axs[1].set_ylabel(f'Dependency', fontdict={'fontsize': LABEL_SIZE}, x=0.05, y=-0.15)

    plt.savefig(os.path.join(FIGURE_DIR, 'Fig_4e_pdac_classical_wnt_scatterplots.pdf'))

# figure 4f

def plot_wnt_densities(wnt_genes_to_regress):
    """
    Plot densities of example Wnt gene dependencies

    Args:
        wnt_genes_to_regress (list[str]): A list of genes to make density plots for
    """
    example_df = load_file('wnt_example_table.csv', local_dir=PROCESSED_DIR)
    
    # create the subplot series
    fig, axs = plt.subplots(4, 1, figsize=(30 * mm, 60 * mm), gridspec_kw={'hspace': 0.05})
    plt.subplots_adjust(left=0.32, bottom=0.18, top=0.98, right=0.92)
    axs = axs.flatten()
    
    # plot the densities and obtain common axis limits
    xlims = np.zeros((len(wnt_genes_to_regress), 2))
    for i, g in enumerate(wnt_genes_to_regress):
        wnt_dependency_subset = example_df[example_df['Dependency'] == g]
        
        sns.kdeplot(
            wnt_dependency_subset,
            x='DependencyGeneEffect',
            hue='SampleSet', palette=sample_set_palette, fill=True,
            linewidth=0.5,
            common_norm=False,
            legend=False,
            ax=axs[i],
        )
        axs[i].axvline(0, linestyle='dotted', color='black', zorder=-1, linewidth=0.5)
        xlims[i, :] = axs[i].get_xlim()

    minx = np.min(xlims) 
    maxx = np.max(xlims)

    # fix the common axis limits
    for i, g in enumerate(wnt_genes_to_regress):
        axs[i].set_xlim(minx, maxx)
        axs[i].set_yticks([])
        axs[i].set_ylabel(g.split(' ')[0], rotation=0, labelpad=15, horizontalalignment='center', fontdict={'size': LABEL_SIZE})

        if i < len(wnt_genes_to_regress) - 1:
            axs[i].set(xticks=[])
            axs[i].set_xlabel('')
            
    axs[3].set_xticks([-3, -2, -1, 0, 1], [-3, -2, -1, 0, 1], fontdict={'size': TICK_SIZE})
    axs[3].set_xlabel('Dependency', fontdict={'size': LABEL_SIZE})
    
    axs[0].text(-0.7, 0.75, 'Organoid', fontdict={'size': ANNOT_SIZE, 'color': sample_set_palette['DepMap Organoid']}, ha='right')
    axs[0].text(0.5, 0.75, '2D', fontdict={'size': ANNOT_SIZE, 'color': sample_set_palette['DepMap 2D']}, ha='left')
    
    plt.savefig(os.path.join(FIGURE_DIR, 'Fig_4f_wnt_dependency_densities.pdf'))
    


def main():
    """
    Generate all plots related to dependencies associated with PDAC-classical expression
    """
    # figure 4d
    plot_pdac_volcano(fdr_threshold=0.05, manual_annotate=MANUAL_ANNOTATE)
    
    # figure 4e
    plot_wnt_regressions(['MESD (23184)', 'WLS (79971)', 'TCF7L2 (6934)', 'FZD5 (7855)'])
    
    # figure 4f
    plot_wnt_densities(['MESD (23184)', 'WLS (79971)', 'TCF7L2 (6934)', 'FZD5 (7855)'])

    
if __name__ == "__main__":
    main()    
