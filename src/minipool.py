import os
import pandas as pd
import numpy as np
from figure_utils import *
from constants import *
from data_utils import *

import seaborn as sns
import matplotlib.pyplot as plt
import upsetplot

mpl.style.use(os.path.join(ASSETS_DIR, 'stylesheet.mplstyle'))

MANUAL_ANNOTATE = True

# figure 5d

def plot_growth_format_volcano(fdr_threshold=0.05, manual_annotate=MANUAL_ANNOTATE):
    """
    Plot differences in dependency across growth formats

    Args:
        fdr_threshold (float, optional): Minimum significance value. Defaults to 0.05.
        manual_annotate (bool, optional): Flag specifying whether to annotate points manually. Defaults to MANUAL_ANNOTATE.
    """
    # load data
    df = load_file('minipool_growth_format_tests.csv', local_dir=PROCESSED_DIR, index_col=0)

    df['-log10(p)'] = -np.log10(df['GrowthPValue'])
    p_boundary_estimate = approximate_fdr_in_p_units(
        df, fdr_threshold=fdr_threshold, pcolumn='GrowthPValue', fdr_column='GrowthFDR'
    )

    plt.figure(figsize=(55 * mm, 40 * mm))
    plt.subplots_adjust(right=0.95, bottom=0.22, top=0.85, left=0.15)

    xlabel='GrowthCoefficient'
    ylabel='-log10(p)'

    # make plot
    sns.scatterplot(
        df[~df['AnyDependency']],
        x=xlabel, y=ylabel,
        s=4,
        hue='RepresentativeGeneset',
        palette=mp_gs_pal,
        alpha=0.3,
        legend=False
    )

    sns.scatterplot(
        df[df['AnyDependency']],
        x=xlabel, y=ylabel,
        s=4,
        hue='RepresentativeGeneset',
        palette=mp_gs_pal,
        alpha=0.8,
        legend=False
    )
    plt.xlabel('$\Delta$ Dependency\n(Dome - Plastic)', fontdict={'size': LABEL_SIZE})
    plt.ylabel('-log10(p)', fontdict={'size': LABEL_SIZE})
    plt.title('Growth Format Differences', fontdict={'size': TITLE_SIZE})
    plt.xticks(plt.gca().get_xticks()[1:-1], labels=plt.gca().get_xticklabels()[1:-1], fontdict={'size': TICK_SIZE})
    plt.yticks(plt.gca().get_yticks()[1:-1], labels=plt.gca().get_yticklabels()[1:-1], fontdict={'size': TICK_SIZE})
    plt.axhline(-np.log10(p_boundary_estimate), linestyle='dashed', alpha=0.5, color='tab:red', linewidth=0.5)
    plt.axvline(0, linestyle='dotted', alpha=0.8, color='black', linewidth=0.5)
    plt.axhline(0, linestyle='dotted', alpha=0.8, color='black', linewidth=0.5)

    # add annotations for specific genes
    annots = df.loc[['PTK2', 'ITGB1', 'ITGA3', 'PRKCA', 'CTNNA1', 'ITGAV', 'NCKAP1'], :]

    if manual_annotate:
        manually_annotate(plt.gca(), "PTK2", annots, 
                          xlabel, ylabel, (0.02, 0.1), ha='left', va='center',
                          arrowprops=dict(arrowstyle='-', color='black', alpha=0, linewidth=0.25, shrinkA=0, shrinkB=1))
        manually_annotate(plt.gca(), "ITGB1", annots, 
                          xlabel, ylabel, (0.02, -0.1), ha='left', va='center',
                          arrowprops=dict(arrowstyle='-', color='black', alpha=0, linewidth=0.25, shrinkA=0, shrinkB=1))
        manually_annotate(plt.gca(), "ITGA3", annots, 
                          xlabel, ylabel, (0.02, 0), ha='left', va='center',
                          arrowprops=dict(arrowstyle='-', color='black', alpha=0, linewidth=0.25, shrinkA=0, shrinkB=1))
        manually_annotate(plt.gca(), "PRKCA", annots, 
                          xlabel, ylabel, (-0.02, 0), ha='right', va='center',
                          arrowprops=dict(arrowstyle='-', color='black', alpha=0, linewidth=0.25, shrinkA=0, shrinkB=0.5))
        manually_annotate(plt.gca(), "CTNNA1", annots, 
                          xlabel, ylabel, (0.02, 0), ha='left', va='center',
                          arrowprops=dict(arrowstyle='-', color='black', alpha=0, linewidth=0.25, shrinkA=0, shrinkB=0))
        manually_annotate(plt.gca(), "ITGAV", annots, 
                          xlabel, ylabel, (-0.02, 0), ha='right', va='center',
                          arrowprops=dict(arrowstyle='-', color='black', alpha=0, linewidth=0.25, shrinkA=0, shrinkB=0))
        manually_annotate(plt.gca(), "NCKAP1", annots, 
                          xlabel, ylabel, (-0.02, 0), ha='right', va='center',
                          arrowprops=dict(arrowstyle='-', color='black', alpha=0, linewidth=0.25, shrinkA=0, shrinkB=0))
    else:
        texts = [plt.gca().text(r[xlabel], r[ylabel], i, ha='center', va='center', fontsize=ANNOT_SIZE) for i,r in annots.iterrows()]
    
    plt.savefig(os.path.join(FIGURE_DIR, 'Fig_5d_minipool_growth_format_volcano.pdf'))
    
# figure 5e

def plot_serum_status_volcano(fdr_threshold=0.05, manual_annotate=MANUAL_ANNOTATE):
    """
    Plot differences in dependency across media conditions

    Args:
        fdr_threshold (float, optional): Minimum significance value. Defaults to 0.05.
        manual_annotate (bool, optional): Flag specifying whether to annotate points manually. Defaults to MANUAL_ANNOTATE.
    """
    # load data
    df = load_file('minipool_serum_status_tests.csv', local_dir=PROCESSED_DIR, index_col=0)

    df['-log10(p)'] = -np.log10(df['MediaPValue'])
    p_boundary_estimate = approximate_fdr_in_p_units(
        df, fdr_threshold=fdr_threshold, pcolumn='MediaPValue', fdr_column='MediaFDR'
    )

    plt.figure(figsize=(55 * mm, 40 * mm))
    plt.subplots_adjust(right=0.95, bottom=0.22, top=0.85, left=0.15)

    xlabel='MediaCoefficient'
    ylabel='-log10(p)'

    # make plot
    sns.scatterplot(
        df[~df['AnyDependency']],
        x=xlabel, y=ylabel,
        s=4,
        hue='RepresentativeGeneset',
        palette=mp_gs_pal,
        alpha=0.3,
        legend=False
    )

    sns.scatterplot(
        df[df['AnyDependency']],
        x=xlabel, y=ylabel,
        s=4,
        hue='RepresentativeGeneset',
        palette=mp_gs_pal,
        alpha=0.8,
        legend=False
    )
    plt.xlabel('$\Delta$ Dependency\n(OPAC - RPMI/FBS)', fontdict={'size': LABEL_SIZE})
    plt.ylabel('-log10(p)', fontdict={'size': LABEL_SIZE})
    plt.title('Culture Medium Differences', fontdict={'size': TITLE_SIZE})
    plt.xticks(plt.gca().get_xticks()[1:-1], labels=plt.gca().get_xticklabels()[1:-1], fontdict={'size': TICK_SIZE})
    plt.yticks(plt.gca().get_yticks()[1:-1], labels=plt.gca().get_yticklabels()[1:-1], fontdict={'size': TICK_SIZE})
    plt.axhline(-np.log10(p_boundary_estimate), linestyle='dashed', alpha=0.5, color='tab:red', linewidth=0.5)
    plt.axvline(0, linestyle='dotted', alpha=0.8, color='black', linewidth=0.5)
    plt.axhline(0, linestyle='dotted', alpha=0.8, color='black', linewidth=0.5)

    # annotate genes
    annots = df.loc[['PTK2', 'RICTOR', 'CTNNA1', 'SPTLC1', 'SPTLC2', 'PPP2R2D', 'PPP1R15B', 'SQLE'], :]

    if manual_annotate:
        manually_annotate(plt.gca(), "PTK2", annots, 
                          xlabel, ylabel, (0.02, 0), ha='left', va='center',
                          arrowprops=dict(arrowstyle='-', color='black', alpha=0, linewidth=0.25, shrinkA=0, shrinkB=1))
        manually_annotate(plt.gca(), "RICTOR", annots, 
                          xlabel, ylabel, (0.02, -0.3), ha='center', va='center',
                          arrowprops=dict(arrowstyle='-', color='black', alpha=0, linewidth=0.25, shrinkA=0, shrinkB=1))
        manually_annotate(plt.gca(), "CTNNA1", annots, 
                          xlabel, ylabel, (0.01, 0.25), ha='left', va='center',
                          arrowprops=dict(arrowstyle='-', color='black', alpha=0, linewidth=0.25, shrinkA=0, shrinkB=1))
        manually_annotate(plt.gca(), "SPTLC2", annots, 
                          xlabel, ylabel, (0, 0.25), ha='left', va='center',
                          arrowprops=dict(arrowstyle='-', color='black', alpha=0, linewidth=0.25, shrinkA=0, shrinkB=0.5))
        manually_annotate(plt.gca(), "SPTLC1", annots, 
                          xlabel, ylabel, (0, 0.6), ha='left', va='top',
                          arrowprops=dict(arrowstyle='-', color='black', alpha=1, linewidth=0.25, shrinkA=8, shrinkB=0))
        manually_annotate(plt.gca(), "SQLE", annots, 
                          xlabel, ylabel, (-0.02, 0), ha='right', va='center',
                          arrowprops=dict(arrowstyle='-', color='black', alpha=0, linewidth=0.25, shrinkA=0, shrinkB=0))
        manually_annotate(plt.gca(), "PPP1R15B", annots, 
                          xlabel, ylabel, (-0.12, 0.3), ha='right', va='center',
                          arrowprops=dict(arrowstyle='-', color='black', alpha=1, linewidth=0.25, shrinkA=0, shrinkB=0))
        manually_annotate(plt.gca(), "PPP2R2D", annots, 
                          xlabel, ylabel, (0.02, 0.4), ha='center', va='bottom',
                          arrowprops=dict(arrowstyle='-', color='black', alpha=1, linewidth=0.25, shrinkA=0, shrinkB=0))
    else:
        texts = [plt.gca().text(r[xlabel], r[ylabel], i, ha='center', va='center', fontsize=ANNOT_SIZE) for i,r in annots.iterrows()]

    plt.savefig(os.path.join(FIGURE_DIR, 'Fig_5e_minipool_serum_status_volcano.pdf'))

# figure 5f

def plot_cross_scatter_plot(manual_annotate=MANUAL_ANNOTATE):
    """
    Plot panels of growth format differences against media condition differences in dependency, highlighting
    genesets for in each panel

    Args:
        manual_annotate (bool, optional): Flag specifying whether to annotate points manually. Defaults to MANUAL_ANNOTATE.
    """
    # load data
    mp_growth_format_tests = load_file('minipool_growth_format_tests.csv', local_dir=PROCESSED_DIR, index_col=0)
    mp_serum_status_tests = load_file('minipool_serum_status_tests.csv', local_dir=PROCESSED_DIR, index_col=0)
    mp_library_composition = load_file('minipool_library_composition.csv', index_col=0)
    comparison_cross_df = mp_growth_format_tests.drop(['RepresentativeGeneset', 'AnyDependency'], axis=1).join(
        mp_serum_status_tests.drop(['RepresentativeGeneset'], axis=1)
    ).join(mp_library_composition)
    
    # decide on genes to annotate per geneset
    genes_to_label = [
        ['ITGAV', 'ITGB1', 'ITGB5', 'ITGA3'], 
        ['PTK2', 'PIK3CA', 'PPP1CB', 'PPP1R12A'],
        ['WLS', 'CTNNB1', 'FZD5', 'TCF7L2', 'APC'],
        ['GPX4', 'SPTLC2', 'SPTLC1', 'SQLE']
    ]
    genesets_to_show = ['Integrin', 'Actin Regulation', 'Wnt Signaling', 'Lipid Metabolism']

    # make plot
    fig, axs = plt.subplots(2, 2, figsize=(60 * mm, 50 * mm), gridspec_kw={'hspace': 0.4, 'wspace': 0.2})
    axs = axs.flatten()
    plt.subplots_adjust(bottom=0.22, left=0.2, right=0.9)
    xlabel='GrowthCoefficient'
    ylabel='MediaCoefficient'

    for i, gs in enumerate(genesets_to_show):
        # genes not in highlighted geneset and shows no dependency in any model
        sns.scatterplot(
            comparison_cross_df[~comparison_cross_df['AnyDependency'] & ~(comparison_cross_df[gs])],
            x=xlabel,
            y=ylabel,
            color='tab:gray',
            alpha=0.3,
            s=4,
            ax=axs[i]
        )

        # genes not in highlighted geneset and shows some dependency in at least one model
        sns.scatterplot(
            comparison_cross_df[comparison_cross_df['AnyDependency'] & ~(comparison_cross_df[gs])],
            x=xlabel,
            y=ylabel,
            color='tab:gray',
            alpha=0.8,
            s=4,
            ax=axs[i]
        )

        # genes not in highlighted geneset but shows some dependency in at least one model
        sns.scatterplot(
            comparison_cross_df[~comparison_cross_df['AnyDependency'] & (comparison_cross_df[gs])],
            x=xlabel,
            y=ylabel,
            color=mp_gs_pal[gs],
            alpha=0.3,
            s=4,
            edgecolor='black',
            ax=axs[i]
        )

        # genes in highlighted geneset and shows some dependency in at least one model
        sns.scatterplot(
            comparison_cross_df[comparison_cross_df['AnyDependency'] & (comparison_cross_df[gs])],
            x=xlabel,
            y=ylabel,
            color=mp_gs_pal[gs],
            alpha=0.8,
            s=4,
            edgecolor='black',
            ax=axs[i]
        )

        axs[i].set_title(gs, fontdict={'fontsize': TITLE_SIZE, 'color': mp_gs_pal[gs]}, pad=3)
        axs[i].axhline(0, c='black', linestyle='dotted', alpha=0.5, linewidth=0.5)
        axs[i].axvline(0, c='black', linestyle='dotted', alpha=0.5, linewidth=0.5)
        axs[i].set_xlabel('')
        axs[i].set_xticks([-0.5, 0, 0.5, 1], [-0.5, 0.0, 0.5, 1], fontdict={'fontsize': TICK_SIZE})
        axs[i].set_yticks([])
        axs[i].set_ylabel('')

        # annotate genes
        annots = comparison_cross_df.loc[genes_to_label[i]]

        if gs == "Integrin":
            if manual_annotate:
                manually_annotate(
                    axs[i], 'ITGAV', annots, xlabel, ylabel, (0.2, 0.3), ha='right', va='center',
                    arrowprops=dict(arrowstyle='-', color='black', alpha=1, linewidth=0.25, shrinkA=0, shrinkB=1)
                )
                manually_annotate(
                    axs[i], 'ITGB5', annots, xlabel, ylabel, (0.09, -0.3), ha='left', va='center',
                    arrowprops=dict(arrowstyle='-', color='black', alpha=1, linewidth=0.25, shrinkA=0, shrinkB=1)
                )
                manually_annotate(
                    axs[i], 'ITGB1', annots, xlabel, ylabel, (0, -0.6), ha='left', va='center',
                    arrowprops=dict(arrowstyle='-', color='black', alpha=1, linewidth=0.25, shrinkA=0, shrinkB=1)
                )
                manually_annotate(
                    axs[i], 'ITGA3', annots, xlabel, ylabel, (-0.12, 0.45), ha='left', va='center',
                    arrowprops=dict(arrowstyle='-', color='black', alpha=1, linewidth=0.25, shrinkA=0, shrinkB=1)
                )
            else:
                texts = [axs[i].text(r[xlabel], r[ylabel], lb, ha='center', va='center', fontsize=ANNOT_SIZE) for lb,r in annots.iterrows()]
        if gs == "Actin Regulation":
            if manual_annotate:
                manually_annotate(
                    axs[i], 'PPP1CB', annots, xlabel, ylabel, (-0.18, -0.08), ha='right', va='center',
                    arrowprops=dict(arrowstyle='-', color='black', alpha=1, linewidth=0.25, shrinkA=0, shrinkB=1)
                )
                manually_annotate(
                    axs[i], 'PPP1R12A', annots, xlabel, ylabel, (-0.02, -0.2), ha='center', va='center',
                    arrowprops=dict(arrowstyle='-', color='black', alpha=1, linewidth=0.25, shrinkA=2.5, shrinkB=1)
                )
                manually_annotate(
                    axs[i], 'PTK2', annots, xlabel, ylabel, (0.2, 0.15), ha='right', va='center',
                    arrowprops=dict(arrowstyle='-', color='black', alpha=0, linewidth=0.25, shrinkA=0, shrinkB=1)
                )
                manually_annotate(
                    axs[i], 'PIK3CA', annots, xlabel, ylabel, (-0.15, 0), ha='right', va='center',
                    arrowprops=dict(arrowstyle='-', color='black', alpha=1, linewidth=0.25, shrinkA=0, shrinkB=1)
                )
            else:
                texts = [axs[i].text(r[xlabel], r[ylabel], lb, ha='center', va='center', fontsize=ANNOT_SIZE) for lb,r in annots.iterrows()]
        if gs == "Wnt Signaling":
            if manual_annotate:
                manually_annotate(
                    axs[i], 'WLS', annots, xlabel, ylabel, (-0.4, -0.5), ha='right', va='center',
                    arrowprops=dict(arrowstyle='-', color='black', alpha=1, linewidth=0.25, shrinkA=0, shrinkB=1)
                )
                manually_annotate(
                    axs[i], 'CTNNB1', annots, xlabel, ylabel, (-0.38, 0.7), ha='center', va='center',
                    arrowprops=dict(arrowstyle='-', color='black', alpha=1, linewidth=0.25, shrinkA=0, shrinkB=1)
                )
                manually_annotate(
                    axs[i], 'FZD5', annots, xlabel, ylabel, (0.25, 0.05), ha='left', va='center',
                    arrowprops=dict(arrowstyle='-', color='black', alpha=1, linewidth=0.25, shrinkA=0, shrinkB=1)
                )
                manually_annotate(
                    axs[i], 'TCF7L2', annots, xlabel, ylabel, (0.2, -0.55), ha='left', va='center',
                    arrowprops=dict(arrowstyle='-', color='black', alpha=1, linewidth=0.25, shrinkA=0, shrinkB=1)
                )
                manually_annotate(
                    axs[i], 'APC', annots, xlabel, ylabel, (0.33, -0.3), ha='left', va='center',
                    arrowprops=dict(arrowstyle='-', color='black', alpha=1, linewidth=0.25, shrinkA=0, shrinkB=1)
                )
            else:
                texts = [axs[i].text(r[xlabel], r[ylabel], lb, ha='center', va='center', fontsize=ANNOT_SIZE) for lb,r in annots.iterrows()]
        if gs == "Lipid Metabolism":
            if manual_annotate:
                manually_annotate(
                    axs[i], 'GPX4', annots, xlabel, ylabel, (0, 0.15), ha='center', va='center',
                    arrowprops=dict(arrowstyle='-', color='black', alpha=0, linewidth=0.25, shrinkA=0, shrinkB=1)
                )
                manually_annotate(
                    axs[i], 'SPTLC2', annots, xlabel, ylabel, (0.2, 0.05), ha='left', va='center',
                    arrowprops=dict(arrowstyle='-', color='black', alpha=1, linewidth=0.25, shrinkA=0, shrinkB=1)
                )
                manually_annotate(
                    axs[i], 'SPTLC1', annots, xlabel, ylabel, (0.28, -0.05), ha='left', va='center',
                    arrowprops=dict(arrowstyle='-', color='black', alpha=1, linewidth=0.25, shrinkA=0, shrinkB=1)
                )
                manually_annotate(
                    axs[i], 'SQLE', annots, xlabel, ylabel, (-0.2, 0), ha='right', va='center',
                    arrowprops=dict(arrowstyle='-', color='black', alpha=1, linewidth=0.25, shrinkA=0, shrinkB=1)
                )
            else:
                texts = [axs[i].text(r[xlabel], r[ylabel], lb, ha='center', va='center', fontsize=ANNOT_SIZE) for lb,r in annots.iterrows()]

    axs[0].set_yticks([-0.5, 0, 0.5], [-0.5, 0.0, 0.5], fontdict={'fontsize': TICK_SIZE})
    axs[2].set_yticks([-0.5, 0, 0.5], [-0.5, 0.0, 0.5], fontdict={'fontsize': TICK_SIZE})
    axs[0].set_xticks([])
    axs[0].set_xlabel('')
    axs[1].set_xticks([])
    axs[1].set_xlabel('')
    axs[2].set_xlabel('')

    axs[0].set_ylabel('$\Delta$ Dependency\n(OPAC - RPMI/FBS)', fontdict={'fontsize': LABEL_SIZE}, y=-0.22)
    axs[2].set_xlabel('$\Delta$ Dependency\n(Dome - Plastic)', fontdict={'size': LABEL_SIZE}, x=1.1)
    plt.savefig(os.path.join(FIGURE_DIR, 'Fig_5e_minipool_cross_scatter_paneled.pdf'))
    
# figure 5g
    
def plot_geneset_variance_scatter(manual_annotate=MANUAL_ANNOTATE):
    """
    Plot F-test results for growth formats against media conditions.

    Args:
        manual_annotate (bool, optional): Flag specifying whether to annotate points manually. Defaults to MANUAL_ANNOTATE.
    """
    mp_geneset_variance_df = load_file('minipool_geneset_variance_tests.csv', local_dir=PROCESSED_DIR)
    
    plt.figure(figsize=(50 * mm, 50 * mm))
    plt.subplots_adjust(bottom=0.22, left=0.28, right=0.95)

    # annotate gene sets
    mp_geneset_variance_df.loc[:, 'HitClass'] = 'Insignificant'
    mp_geneset_variance_df.loc[mp_geneset_variance_df['GrowthFDR'] < 0.05, 'HitClass'] = 'Growth'
    mp_geneset_variance_df.loc[mp_geneset_variance_df['MediaFDR'] < 0.05, 'HitClass'] = 'Media'

    xlabel = 'GrowthVariance'
    ylabel = 'MediaVariance'

    # make plot
    sns.scatterplot(
        mp_geneset_variance_df,
        x=xlabel, y=ylabel,
        hue='HitClass', palette={'Insignificant': 'tab:gray', 'Growth': 'tab:blue', 'Media': 'tab:red'},
        s=8,
        legend=False
    )

    plt.ylabel('Variance[$\Delta$ Dependency]\nCulture Medium')
    plt.xlabel('Variance[$\Delta$ Dependency]\nGrowth Format')

    minx, maxx = plt.xlim()
    miny, maxy = plt.ylim()
    plt.xlim(0, maxx)
    plt.ylim(0, maxy)

    # annotate significant points
    annots = pd.concat([
        mp_geneset_variance_df[(mp_geneset_variance_df['GrowthFDR'] < 0.05) | (mp_geneset_variance_df['MediaFDR'] < 0.05)],
        mp_geneset_variance_df[mp_geneset_variance_df['Geneset'].isin(['Wnt Signaling', 'All test genes'])]
    ], axis=0).drop_duplicates().set_index('Geneset')

    if manual_annotate:
        manually_annotate(
            plt.gca(), 'Lipid Metabolism', annots, xlabel, ylabel, (0, -0.005), ha='center', va='center',
            arrowprops=dict(arrowstyle='-', color='black', alpha=0, linewidth=0.25, shrinkA=0, shrinkB=1)
        )
        manually_annotate(
            plt.gca(), 'Actin Regulation', annots, xlabel, ylabel, (0, 0.01), ha='center', va='center',
            arrowprops=dict(arrowstyle='-', color='black', alpha=1, linewidth=0.25, shrinkA=0, shrinkB=1)
        )
        manually_annotate(
            plt.gca(), 'Integrin', annots, xlabel, ylabel, (-0.005, 0), ha='right', va='center',
            arrowprops=dict(arrowstyle='-', color='black', alpha=1, linewidth=0.25, shrinkA=0, shrinkB=1)
        )
        manually_annotate(
            plt.gca(), 'All test genes', annots, xlabel, ylabel, (0.007, 0), ha='left', va='center',
            arrowprops=dict(arrowstyle='-', color='black', alpha=1, linewidth=0.25, shrinkA=0, shrinkB=1)
        )
        manually_annotate(
            plt.gca(), 'Wnt Signaling', annots, xlabel, ylabel, (0, -0.01), ha='center', va='center',
            arrowprops=dict(arrowstyle='-', color='black', alpha=1, linewidth=0.25, shrinkA=0, shrinkB=1)
        )
    else:
        texts = [plt.gca().text(r[xlabel], r[ylabel], i, ha='center', va='center', fontsize=ANNOT_SIZE) for i,r in annots.loc[['Lipid Metabolism', 'Actin Regulation', 'Integrin', 'All test genes', 'Wnt Signaling'], :].iterrows()]

    legend = plt.legend(
        handles=[
            Line2D([], [], color='tab:gray', marker = 'o', linestyle='None', markersize=4, markeredgewidth=0, label = 'Insigificant'),
            Line2D([], [], color='tab:blue', marker = 'o', linestyle='None', markersize=4, markeredgewidth=0, label = 'High variance\nacross format'),
            Line2D([], [], color='tab:red', marker = 'o', linestyle='None', markersize=4, markeredgewidth=0, label = 'High variance\nacross medium')
        ], handletextpad=0.2, prop={'size': ANNOT_SIZE}
    )
    plt.setp(legend.get_title(), fontsize=ANNOT_SIZE)

    plt.savefig(os.path.join(FIGURE_DIR, 'Fig_5g_minipool_geneset_variance_scatter.pdf'))
    
# figure ED11a
    
def plot_minipool_library_composition():
    """
    Make an upset plot showing the membership of genes in the minipool with selected gene sets
    """
    test_gene_annotations = load_file('minipool_library_composition.csv', index_col=0)
    fig = plt.figure(figsize=(180 * mm, 90 * mm))
    mtx_up = test_gene_annotations.drop(['Other', 'All test genes', 'Actin Nucleation', 'Gap Junction'], axis=1, errors='ignore')
    upsetplot.plot(
        upsetplot.from_indicators(mtx_up), 
        subset_size='count', sort_categories_by='-cardinality', sort_by='cardinality',
        fig=fig, element_size=10
    )
    plt.savefig(os.path.join(FIGURE_DIR, 'Fig_ED11a_minipool_library_upsetplot.pdf'), bbox_inches='tight')

# figure ED11b

def plot_minipool_screen_qc():
    """
    Plot quality control metrics for the minipool screens
    """
    # load QC data
    mp_screen_metadata = load_file('minipool_screen_metadata.csv')
    
    # group screens into colors by their model ID
    cl_qc_colormap = dict()
    cl_colors = sns.color_palette('tab20', n_colors=20)
    for i, dmid in enumerate(mp_screen_metadata['DepMap ID'].unique()):
        cl_qc_colormap[dmid] = cl_colors[i*2]
        cl_qc_colormap[dmid + '_'] = cl_colors[(i*2) + 1]

    cl_qc_patches = [Patch(color=cl_qc_colormap[dmid], label=dmid) for dmid in mp_screen_metadata['DepMap ID'].unique()]

    fig, axs = plt.subplots(1, 3, figsize=(120 * mm, 80 * mm), sharey=True)

    # plot total recovered reads
    sns.barplot(
        mp_screen_metadata,
        y='cell_line_name',
        x='Total reads',
        hue='DepMap ID',
        palette=cl_qc_colormap,
        dodge=False,
        ax=axs[0]
    )
    axs[0].set_ylabel('Screen Name')
    axs[0].set_xlabel('Total Reads')
    axs[0].axvline(5e6, color='black', linestyle='dashed', linewidth=0.5)
    axs[0].legend_.remove()
    axs[0].set_title('Total Readcounts')

    # plot NNMD scores
    sns.barplot(
        mp_screen_metadata,
        y='cell_line_name',
        x='NNMD',
        hue='DepMap ID',
        palette=cl_qc_colormap,
        dodge=False,
        ax=axs[1],
    )
    axs[1].set_xlabel('NNMD')
    axs[1].axvline(-1.5, color='black', linestyle='dashed', linewidth=0.5)
    axs[1].legend_.remove()
    axs[1].set_title('Screen NNMD')

    # pot replicate correlation
    sns.barplot(
        mp_screen_metadata,
        y='cell_line_name',
        x='Replicate LFC Correlation',
        hue='DepMap ID',
        palette=cl_qc_colormap,
        dodge=False,
        ax=axs[2]
    )
    axs[2].set_xlabel("Pearson's r")
    axs[2].axvline(0.2, color='black', linestyle='dashed', linewidth=0.5)
    axs[2].legend_.remove()
    axs[2].set_title('Replicate LFC Correlation')

    plt.legend(title='DepMap ID', handles=cl_qc_patches, loc='lower left', bbox_to_anchor=(1.04, 0.5))
            
    plt.savefig(os.path.join(FIGURE_DIR, 'Fig_ED11b_minipool_screen_qc.pdf'), bbox_inches='tight')
    
# figure ED11c

def plot_minipool_genomewide_comparison():
    """
    Plot gene effects in the genome-wide screens against their minipool screens corresponding screening conditions.
    """
    # load data
    mp_sample_metadata = load_file('minipool_sequence_metadata.csv')
    mp_screen_metadata = load_file('minipool_screen_metadata.csv')
    mp_chronos_gene_effect = load_file('minipool_chronos_gene_effects.csv', local_dir=PROCESSED_DIR, index_col=0)
    mp_guide_metadata = load_file('minipool_guide_metadata.csv')
    screen_metadata = load_file('screen_metadata.csv', index_col=0)
    screen_gene_effect = load_file('screen_gene_effect.csv', index_col=0)
    
    achid_to_native_growth = mp_sample_metadata.query('`Media condition` == "Native"').drop_duplicates(
        ['DepMap ID', 'Original Growth Pattern']
    ).set_index('DepMap ID')['Original Growth Pattern'].to_dict()
    
    # map native screen conditions for each model
    native_screens = mp_screen_metadata[
        (mp_screen_metadata['Media condition'] == "Native") &
        (((mp_screen_metadata['Original Growth Pattern'] == "Adherent") & 
          (mp_screen_metadata['Growth Format'] == "Plastic")) |
         ((mp_screen_metadata['Original Growth Pattern'] == "Organoid") & 
          (mp_screen_metadata['Growth Format'] == "Dome")))
    ].set_index('cell_line_name')['DepMap ID']

    screen_id_to_common_depmap_id = native_screens.reset_index().merge(
        screen_metadata.reset_index().loc[:, ['ModelID', 'ScreenID']], 
        left_on='DepMap ID', right_on='ModelID'
    ).set_index('ScreenID')['ModelID']
    
    # subset to only the genes in the minipool library
    common_minipool_gene_effect = mp_chronos_gene_effect.reindex(
        index=native_screens.index.tolist()
    ).rename(native_screens).dropna()

    common_genomewide_gene_effect = screen_gene_effect.groupby(screen_id_to_common_depmap_id).median()
    common_genomewide_gene_effect = common_genomewide_gene_effect.rename(
        {c:c.split(' ')[0] for c in common_genomewide_gene_effect.columns}, axis=1
    ).reindex_like(common_minipool_gene_effect).dropna(how = "all")
    
    # annotate positive and negative control genes
    pos_ctrl_genes = mp_guide_metadata.query('Gene_category == "Positive_control"')['Gene'].unique().tolist()
    neg_ctrl_genes = mp_guide_metadata.query('Gene_category == "Negative_control"')['Gene'].unique().tolist()
    test_genes = mp_guide_metadata.query('Gene_category == "Test_gene"')['Gene'].unique().tolist()

    gene_class = pd.concat([
        pd.Series('Nonessential', index=neg_ctrl_genes),
        pd.Series('Essential', index=pos_ctrl_genes),
        pd.Series('Test', index=test_genes)
    ]).rename('Gene Category')
    
    # plot panels of scatterplots
    ncols = 2
    nrows = common_minipool_gene_effect.shape[0] // ncols
    if common_minipool_gene_effect.shape[0] % ncols != 0:
        nrows += 1
    fig, axs = plt.subplots(nrows, ncols, figsize=((ncols * 40 * mm) + (35 * mm), nrows * 40 * mm), sharey=True)
    axs=axs.flatten()

    gene_palette = {'Nonessential': '#0343DF', 'Essential': '#C20078', 'Test': 'tab:gray'}

    model_comparison_table = []
    for i, cl in enumerate(common_minipool_gene_effect.index.tolist()):
        if cl not in common_genomewide_gene_effect.index:
            break

        print(i, cl)
        model_comparison_df = pd.concat([
            pd.Series(cl, index=gene_class.index.tolist(), name='DepMap ID'),
            common_genomewide_gene_effect.loc[cl, gene_class.index.tolist()].rename(f'Genome-wide'),
            common_minipool_gene_effect.loc[cl, gene_class.index.tolist()].rename(f'Minipool'),
            gene_class
        ], axis=1).dropna().reset_index().rename({'index': 'Gene'}, axis=1)
        
        sns.scatterplot(
            model_comparison_df,
            x=f'Genome-wide',
            y=f'Minipool',
            hue='Gene Category', palette=gene_palette,
            s=4,
            ax=axs[i], legend=False
        )

        axs[i].axhline(0, linestyle='dotted', color='black', linewidth=0.5)
        axs[i].axvline(0, linestyle='dotted', color='black', linewidth=0.5)
        axs[i].axline((0, 0), slope=1, linestyle='dotted', color='black', linewidth=0.5)
        axs[i].set_ylabel(f'Minipool')
        axs[i].set_xlabel(f'Genome-wide')
        axs[i].set_title(f'{cl}')
        axs[i].set_xlim(-5, 1)

        # annotate comparisons with the correlations across screens
        axs[i].text(
            0.05, 0.91, 'r = {:.3f}'.format(scipy.stats.pearsonr(model_comparison_df[f'Genome-wide'], model_comparison_df[f'Minipool'])[0]), ha='left', transform=axs[i].transAxes, fontdict={'size': LABEL_SIZE}
        )
        model_comparison_table.append(model_comparison_df)
    model_comparison_table = pd.concat(model_comparison_table, axis=0)

    handles = []
    for gc in gene_palette:
        handle = Line2D([], [], color=gene_palette[gc], marker = 'o', linestyle='None', markersize=4, markeredgewidth=0, label = gc)
        handles.append(handle)    
    plt.legend(handles = handles, loc='lower left', bbox_to_anchor=(1.05, 0))
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURE_DIR, 'minipool_vs_genomewide_scatterplots.pdf'), bbox_inches='tight')
    
    model_comparison_table.rename({
        'Genome-wide': 'Genome-wide Dependency', 'Minipool': 'Minipool Dependency'
    }, axis=1).loc[:, [
        'DepMap ID', 'Gene', 'Genome-wide Dependency', 'Minipool Dependency', 'Gene Category'
    ]].to_csv(os.path.join(PROCESSED_DIR, 'minipool_vs_genomewide_screen_table.csv'), index=False)

# figures ED12ab

def get_pairs(split_field, metadata, 
              matched_fields = ['DepMap ID', 'Original Growth Pattern', 'Media condition', 'Serum status', 'Growth Format']):
    """
    Identify all model pairs that are distinguished by a change in only growth format or only media condition.

    Args:
        split_field (str): A column label which will have differential values between model pairs
        metadata (pandas.DataFrame): A table of screen metadata
        matched_fields (list, optional): A list of column attributes that must be identical between pairs (except `split_field`, which will be ignored). Defaults to ['DepMap ID', 'Original Growth Pattern', 'Media condition', 'Serum status', 'Growth Format'].

    Returns:
        pandas.DataFrame: A table of screen pairs
    """
    matched_fields = list(set(matched_fields) - {split_field})
    
    pairs = dict()
    split_field_options = sorted(list(metadata[split_field].unique()))
    
    # indices will be split_field_options[0], values will be split_field_options[1]
    for i, r in metadata.iterrows():
        if i not in pairs.keys() and i not in pairs.values():
            search = metadata.query(
                ' & '.join([f'(`{m}` == "{r[m]}")' for m in matched_fields]) + f' & (`{split_field}` != "{r[split_field]}")'
            )
            if search.shape[0] == 1:
                if r[split_field] == split_field_options[0]:
                    pairs[i] = search.index[0]
                elif r[split_field] == split_field_options[1]:
                    pairs[search.index[0]] = i
    pairs = pd.Series(pairs).rename(split_field_options[1])
    pairs.index.name = split_field_options[0]
    return pairs.reset_index()

def plot_paired_ttest(genes, gene_effect_matrix, grouping_df, figsize=None, p_annotation=None, 
                      rotate_group_label=0, metadata=None, hue=None, palette='hls', xticklabelpad=-0.05):
    """
    Plot paired dependency data for select genes across screen pairs

    Args:
        genes (list[str]): The list of genes to show
        gene_effect_matrix (pandas.DataFrame): A matrix (screen x gene) of gene effect scores
        grouping_df (pandas.DataFrame): A table of screen pairs, output by `get_pairs()`

    Returns:
        pandas.DataFrame: A table of all values underlying the paired boxplots
    """
    assert grouping_df.shape[1] == 2, 'grouping_df must have exactly 2 columns'
    
    # subset the gene effect matrix by groups
    group1_name = grouping_df.columns[0]
    group2_name = grouping_df.columns[1]
    
    group1_mtx = gene_effect_matrix.loc[grouping_df[group1_name].tolist(), genes]
    group2_mtx = gene_effect_matrix.loc[grouping_df[group2_name].tolist(), genes]
    
    # join values together into a longform matrix
    long_df = pd.concat([
        group1_mtx.melt(ignore_index=False).reset_index().assign(Group=group1_name),
        group2_mtx.melt(ignore_index=False).reset_index().assign(Group=group2_name)
    ], axis=0).rename({'index': 'cell_line_name', 'variable': 'gene'}, axis=1)
    long_df['gene_group'] = long_df['gene'] + ':' + long_df['Group']
    long_df = long_df.sort_values('gene_group')
    
    # group values by gene knockout and set coordinates that will keep pairs together
    discretized_gene_groups = dict()
    discretized_genes = dict()
    for i, g in enumerate(genes):
        discretized_gene_groups[g + ':' + group1_name] = i*2
        discretized_gene_groups[g + ':' + group2_name] = (i*2) + 1
        discretized_genes[g] = (i*2) + 0.5

    position_to_gene_group = {v:k for k,v in discretized_gene_groups.items()}
    position_to_gene = {v:k for k,v in discretized_genes.items()}
    
    long_df['discretized_x'] = long_df['gene_group'].replace(discretized_gene_groups)
    if metadata is not None:
        long_df = long_df.merge(metadata, left_on='cell_line_name', right_index=True)
        
    if figsize is None:
        figsize = (6, 4)
    plt.figure(figsize=figsize)
    
    # make background boxplots
    ax = sns.boxplot(
        long_df,
        x='discretized_x', y='value', showfliers=False, 
        medianprops={'color': 'red'},
        boxprops={"facecolor": (1, 1, 1, 0.1), 'edgecolor': (0, 0, 0, 0.5)},
        whiskerprops={'color': (0, 0, 0, 0.5)},
        capprops={'color': (0, 0, 0, 0.5)}
    )
    
    # plot points on top of boxplots
    if metadata is not None and hue is not None:
        if type(palette) == dict:
            element_to_color = palette
        else:
            palette = sns.color_palette(palette, len(long_df[hue].unique()))
            element_to_color = {el: palette[i] for i, el in enumerate(long_df[hue].unique())}
        ax = sns.scatterplot(
            long_df,
            x='discretized_x', y='value', hue=hue, palette=element_to_color
        )
    else:
        ax = sns.scatterplot(
            long_df,
            x='discretized_x', y='value', color='black'
        )
    
    # join paired data with lines
    for g in genes:
        for i, pair in grouping_df.iterrows():
            if metadata is not None and hue is not None:
                color = element_to_color[metadata.loc[pair[group1_name], hue]]
            else:
                color = 'black'
            
            plt.plot(
                [discretized_gene_groups[g + ':' + group1_name], discretized_gene_groups[g + ':' + group2_name]],
                [gene_effect_matrix.loc[pair[group1_name], g], gene_effect_matrix.loc[pair[group2_name], g]],
                c=color
            )
    
    x_tick_positions = sorted(list(discretized_genes.values()) + list(discretized_gene_groups.values()))
    ax.set_xticks(x_tick_positions)
    
    # annotate gene names along the x-axis with signification values
    if p_annotation is None:
        x_tick_labels = [position_to_gene[x] if x in position_to_gene else position_to_gene_group[x].split(':')[-1]
                         for x in sorted(x_tick_positions)]
    elif type(p_annotation) == pd.Series:
        p_annotation = p_annotation.to_dict()
        x_tick_labels = [position_to_gene[x] + f'\np={p_annotation[position_to_gene[x]]:.3f}' if x in position_to_gene 
                         else position_to_gene_group[x].split(':')[-1]
                         for x in sorted(x_tick_positions)]
        
    ax.set_xticklabels(x_tick_labels)
    for t, position in zip(ax.get_xticklabels(), x_tick_positions):
        if int(position) != position:
            t.set_y(-0.05 + xticklabelpad)
        else:
            t.set_rotation(rotate_group_label)
    ax.tick_params(axis='x', which='major', bottom=False)
    ax.set_xlabel('')
    ax.set_ylabel('Gene Effect')
    
    return long_df

# figure ED12a

def plot_growth_paired_t_boxplots(genelist):
    """
    Plot paired t-test results for growth format comparisons

    Args:
        genelist (list[str]): A list of genes to display
    """
    # load data
    mp_chronos_gene_effect = load_file('minipool_chronos_gene_effects.csv', local_dir=PROCESSED_DIR, index_col=0)
    mp_screen_metadata = load_file('minipool_screen_metadata.csv')
    growth_format_tests = load_file('minipool_growth_format_tests.csv', local_dir=PROCESSED_DIR, index_col=0)
    
    # plot paired boxplots
    growth_format_pairing = get_pairs('Growth Format', metadata=mp_screen_metadata.set_index('cell_line_name'))
    growth_long_df = plot_paired_ttest(
        genelist, 
        mp_chronos_gene_effect, growth_format_pairing, 
        p_annotation=growth_format_tests['GrowthFDR'],
        metadata=mp_screen_metadata.set_index('cell_line_name'), 
        hue='Culture medium', palette={'RPMI/FBS': 'tab:blue', 'OPAC': 'tab:orange'}, 
        figsize=(200 * mm, 60 * mm), 
        rotate_group_label=90, xticklabelpad=-0.12
    )
    plt.title('Chronos Gene Effects: Dome vs Plastic')
    plt.axhline(0, linestyle='dotted', c='black')
    
    plt.savefig(os.path.join(FIGURE_DIR, 'minipool_growth_format_pairs_boxplots.pdf'), bbox_inches='tight')
    
    growth_long_df.rename({'gene': 'Gene', 'value': 'DependencyGeneEffect'}, axis=1).loc[
        :, ['cell_line_name', 'Gene', 'DependencyGeneEffect', 'Growth Format', 'Culture medium']
    ].to_csv(os.path.join(PROCESSED_DIR, 'minipool_growth_format_paired_t_table.csv'), index=False)

# figure ED12b

def plot_serum_paired_t_boxplots(genelist):
    """
    Plot paired t-test results for media comparisons

    Args:
        genelist (list[str]): A list of genes to display
    """
    # load data
    mp_chronos_gene_effect = load_file('minipool_chronos_gene_effects.csv', local_dir=PROCESSED_DIR, index_col=0)
    mp_screen_metadata = load_file('minipool_screen_metadata.csv')
    serum_status_tests = load_file('minipool_serum_status_tests.csv', local_dir=PROCESSED_DIR, index_col=0)
    
    # plot paired boxplots
    serum_status_pairing = get_pairs('Culture medium', metadata=mp_screen_metadata.set_index('cell_line_name'),
                                     matched_fields=['DepMap ID', 'Original Growth Pattern', 'Growth Format'])
    serum_long_df = plot_paired_ttest(
        genelist, 
        mp_chronos_gene_effect, serum_status_pairing, 
        p_annotation=serum_status_tests['MediaFDR'],
        metadata=mp_screen_metadata.set_index('cell_line_name'),
        hue='Growth Format', palette={'Plastic': 'tab:blue', 'Dome': 'tab:orange'},
        figsize=(200 * mm, 60 * mm),
        rotate_group_label=90, xticklabelpad=-0.18
    )
    
    plt.title('Chronos Gene Effects: OPAC vs. RPMI/FBS')
    plt.axhline(0, linestyle='dotted', c='black')
    
    plt.savefig(os.path.join(FIGURE_DIR, 'minipool_serum_status_pairs_boxplots.pdf'), bbox_inches='tight')
    
    serum_long_df.rename({'gene': 'Gene', 'value': 'DependencyGeneEffect'}, axis=1).loc[
        :, ['cell_line_name', 'Gene', 'DependencyGeneEffect', 'Growth Format', 'Culture medium']
    ].to_csv(os.path.join(PROCESSED_DIR, 'minipool_serum_status_paired_t_table.csv'), index=False)

# figure ED12c

def plot_minipool_geneset_boxplots_over_growth():
    """
    Plot genesets and their member genes as boxplots of dependency differences over growth formats
    """
    # load data
    growth_format_tests = load_file('minipool_growth_format_tests.csv', local_dir=PROCESSED_DIR, index_col=0)
    mp_geneset_variance_df = load_file('minipool_geneset_variance_tests.csv', local_dir=PROCESSED_DIR)
    test_gene_annotations = load_file('minipool_library_composition.csv', index_col=0)
    
    # extract values per gene and the genesets they belong to
    df = growth_format_tests.join(
        test_gene_annotations.drop(
            ['Other', 'All test genes', 'Actin Nucleation', 'Gap Junction'], axis=1, errors='ignore'
        ).melt(ignore_index=False).query('value')
    ).reset_index()
    gs_order = mp_geneset_variance_df.sort_values('GrowthVariance', ascending=False)['Geneset'].loc[
        lambda x: x.isin(df['variable'].tolist())
    ].tolist()

    plt.figure(figsize=(70 * mm, 50 * mm))
    
    # recolor non-highlighted genesets
    extended_mp_gs_pal = mp_gs_pal.copy()
    for gs in gs_order:
        if gs not in extended_mp_gs_pal:
            extended_mp_gs_pal[gs] = mp_gs_pal['Other']

    # plot boxplots and individual gene points
    sns.boxplot(
        df, x='GrowthCoefficient', y='variable', 
        order=gs_order,
        showfliers=False, 
        hue='variable', palette=extended_mp_gs_pal, 
        medianprops={'color': 'black'},
        boxprops={"alpha": 0.5},
        whiskerprops={'alpha': 0.5, 'linewidth': 0.5},
        capprops={'alpha': 0.5, 'linewidth': 0.5},
        legend=False
    )
    sns.stripplot(
        df, x='GrowthCoefficient', y='variable', 
        order=gs_order, 
        hue='variable', palette=extended_mp_gs_pal,
        edgecolor='black', linewidth=0.25,
        size=2,
        legend=False,
    )

    plt.xlabel('$\Delta$ Dependency\n(Dome - Plastic)', fontdict={'size': LABEL_SIZE})
    plt.ylabel('Geneset', fontdict={'size': LABEL_SIZE})
    # plt.title('Genesets over growth format', fontdict={'size': TITLE_SIZE})
    plt.axvline(0, linestyle='dotted', color='black', linewidth=0.5, zorder=-10)
    plt.savefig(os.path.join(FIGURE_DIR, 'minipool_growth_format_geneset_boxplots.pdf'), bbox_inches='tight')
    
    df.rename({'variable': 'Geneset'}, axis=1).loc[
        :, ['Gene', 'GrowthCoefficient', 'Geneset']
    ].to_csv(os.path.join(PROCESSED_DIR, 'minipool_growth_format_geneset_boxplot_table.csv'), index=False)

# figure ED12d

def plot_geneset_boxplots_over_serum():
    """
    Plot genesets and their member genes as boxplots of dependency differences over media conditions
    """
    # load data
    serum_status_tests = load_file('minipool_serum_status_tests.csv', local_dir=PROCESSED_DIR, index_col=0)
    mp_geneset_variance_df = load_file('minipool_geneset_variance_tests.csv', local_dir=PROCESSED_DIR)
    test_gene_annotations = load_file('minipool_library_composition.csv', index_col=0)
    
    # extract values per gene and the genesets they belong to
    df = serum_status_tests.join(
        test_gene_annotations.drop(
            ['Other', 'All test genes', 'Actin Nucleation', 'Gap Junction'], axis=1, errors='ignore'
        ).melt(ignore_index=False).query('value')
    ).reset_index()
    gs_order = mp_geneset_variance_df.sort_values('MediaVariance', ascending=False)['Geneset'].loc[
        lambda x: x.isin(df['variable'].tolist())
    ].tolist()

    plt.figure(figsize=(70 * mm, 50 * mm))
    
    # recolor non-highlighted genesets
    extended_mp_gs_pal = mp_gs_pal.copy()
    for gs in gs_order:
        if gs not in extended_mp_gs_pal:
            extended_mp_gs_pal[gs] = mp_gs_pal['Other']

    # plot boxplots and individual gene points
    sns.boxplot(
        df, x='MediaCoefficient', y='variable', 
        order=gs_order,
        showfliers=False, 
        hue='variable', palette=extended_mp_gs_pal, 
        medianprops={'color': 'black'},
        boxprops={"alpha": 0.5},
        whiskerprops={'alpha': 0.5, 'linewidth': 0.5},
        capprops={'alpha': 0.5, 'linewidth': 0.5},
        legend=False
    )
    sns.stripplot(
        df, x='MediaCoefficient', y='variable', 
        order=gs_order, 
        hue='variable', palette=extended_mp_gs_pal,
        edgecolor='black', linewidth=0.25,
        size=2,
        legend=False,
    )

    plt.xlabel('$\Delta$ Dependency\n(OPAC - RPMI/FBS)', fontdict={'size': LABEL_SIZE})
    plt.ylabel('Geneset', fontdict={'size': LABEL_SIZE})
    plt.title('Genesets over culture medium', fontdict={'size': TITLE_SIZE})
    plt.axvline(0, linestyle='dotted', color='black', linewidth=0.5, zorder=-10)
    plt.savefig(os.path.join(FIGURE_DIR, 'minipool_serum_status_geneset_boxplots.pdf'), bbox_inches='tight')
    
    df.rename({'variable': 'Geneset'}, axis=1).loc[
        :, ['Gene', 'MediaCoefficient', 'Geneset']
    ].to_csv(os.path.join(PROCESSED_DIR, 'minipool_serum_status_geneset_boxplot_table.csv'), index=False)
    
def main():
    """
    Generate all plots related to minipool screens
    """
    # # figure 5d
    plot_growth_format_volcano(fdr_threshold=0.05, manual_annotate=MANUAL_ANNOTATE)
    
    # # figure 5e
    plot_serum_status_volcano(fdr_threshold=0.05, manual_annotate=MANUAL_ANNOTATE)
    
    # # figure 5f
    plot_cross_scatter_plot(manual_annotate=MANUAL_ANNOTATE)
    
    # # figure 5g
    plot_geneset_variance_scatter(manual_annotate=MANUAL_ANNOTATE)
    
    # # figure ED11a
    plot_minipool_library_composition()
    
    # # figure ED11b
    plot_minipool_screen_qc()
    
    # figure ED11c
    plot_minipool_genomewide_comparison()
    
    genelist = ['ITGB1', 'PTK2', 'ITGA3', 'PRKCA', 'CTNNA1', 'NCKAP1', 'ITGAV', 'SPTLC2', 'SQLE', 'SPTLC1', 'PPP1R15B', 'PPP2R2D', 'RICTOR', 'CDK4', 'GPX4']
    
    # # figure ED12a
    plot_growth_paired_t_boxplots(genelist)
    
    # # figure ED12b
    plot_serum_paired_t_boxplots(genelist)
    
    # # figure ED12c
    plot_minipool_geneset_boxplots_over_growth()
    
    # # figure ED12d
    plot_geneset_boxplots_over_serum()

    
if __name__ == "__main__":
    main()    
    
