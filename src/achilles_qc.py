import os
import pandas as pd
import numpy as np
from figure_utils import *
from constants import *
from data_utils import *

import seaborn as sns
import matplotlib.pyplot as plt

mpl.style.use(os.path.join(ASSETS_DIR, 'stylesheet.mplstyle'))

def load_data():
    """
    Load all data necessary to make plots about CRISPR screen quality control
    """
    screen_qc_table = load_file('crispr_screen_qc.csv')
    ngs = load_file('crispr_naive_gene_score.csv', index_col=0)
    controls = load_file('crispr_control_genes.csv')
    ness = controls[controls['Category'] == "Nonessential"]['Gene']
    ess = controls[controls['Category'] == "CommonEssential"]['Gene']
    screen_metadata = load_file('screen_metadata.csv', index_col=0)
    screen_gene_effect = load_file('screen_gene_effect.csv', index_col=0)
    screen_dict = screen_qc_table.groupby('Group').apply(lambda x: list(x['ScreenID'])).to_dict()
    screen_dict = {grp: pd.Series(scr) for grp, scr in screen_dict.items()}
    
    return ngs, ness, ess, screen_dict, screen_qc_table, screen_metadata, screen_gene_effect

# figure 1g

def plot_control_separations(ngs, screen_dict, ness, ess):
    """
    Plots the distribution of log-fold changes in reads per gene for non-essential and pan-essential genes across
    different sets of screens.

    Args:
        ngs (pandas.DataFrame): A matrix (screen x gene) of "naive" gene effects, which are average log foldchanges in counts across guides
        screen_dict (dict[str -> pandas.DataFrame]): A mapping of screen subsets to plot
        ness (list[str]): A list of non-essential genes
        ess (list[str]): A list of pan-essential genes
    """
    # create a long-form table of gene effect scores and categorize genes by their control set
    ngs_controls_longform = dict()
    for grp in screen_dict:
        ngs_controls_longform[grp] = pd.concat([
            ngs.reindex(index=screen_dict[grp].values.tolist(), columns=ness).mean().rename('MeanNaiveGE').reset_index().assign(ControlCategory='Nonessential'),
            ngs.reindex(index=screen_dict[grp].values.tolist(), columns=ess).mean().rename('MeanNaiveGE').reset_index().assign(ControlCategory='Essential')
        ])
    
    qc_categories = list(ngs_controls_longform.keys())
    ncat = len(qc_categories)
    fig, axs = plt.subplots(ncat, 1, figsize=(40*mm, 70*mm), sharey=False)
    axs = axs.flatten()
    plt.subplots_adjust(left=0.08, right=0.94, top=0.95, bottom=0.15)

    # plot the gene effect distributions as densities
    xlims = np.zeros((ncat, 2))
    for i, ms in enumerate(qc_categories):
        subset = ngs_controls_longform[ms].reset_index(drop=True)
        sns.kdeplot(subset, x='MeanNaiveGE', hue='ControlCategory', 
                    palette={'Nonessential': '#0343DF', 'Essential': '#C20078'}, fill=True, ax=axs[i])
        axs[i].axhline(0, c='black')
        xlims[i, :] = axs[i].get_xlim()

    minx = -5.25#np.min(xlims) 
    maxx = 1#np.max(xlims)

    # reformat all plots and remove extra legends
    for i, ms in enumerate(qc_categories):
        axs[i].set_xlim(minx, maxx)

        axs[i].spines[['bottom', 'left', 'right', 'top']].set_visible(False)
        axs[i].set(yticks=[])
        axs[i].set_ylabel('')
        axs[i].get_legend().remove()
        if i < ncat-1:
            # axs[i].get_legend().remove()
            axs[i].set(xticks=[])
            axs[i].set_xlabel('')
        else:
            axs[i].set_xlabel('Mean Naive Gene Effect', fontdict={'size': LABEL_SIZE})

    # draw the NNMD components as annotations
    xlim_width = maxx - minx

    for i, ms in enumerate(qc_categories):
        subset = ngs_controls_longform[ms].reset_index(drop=True)
        medians = subset.loc[:, ['ControlCategory', 'MeanNaiveGE']].groupby('ControlCategory').median()
        med_ess = medians.loc['Essential', 'MeanNaiveGE']
        med_ness = medians.loc['Nonessential', 'MeanNaiveGE']
        mad = (subset.query('ControlCategory == "Nonessential"')['MeanNaiveGE'] - 
               subset.query('ControlCategory == "Nonessential"')['MeanNaiveGE'].median()).abs().median()
        mad_minlim = med_ness - mad
        mad_maxlim = med_ness + mad

        med_ess_transform = (med_ess - minx) / xlim_width
        med_ness_transform = (med_ness - minx) / xlim_width
        mad_min_transform = (mad_minlim - minx) / xlim_width
        mad_max_transform = (mad_maxlim - minx) / xlim_width

        # ylabel
        axs[i].text(0, 0.4, ms, horizontalalignment='left', transform=axs[i].transAxes, fontdict={'color': 'black', 'size': LABEL_SIZE})

        # median bars
        axs[i].plot([med_ess, med_ess], [0, 0.2], c='#F5F5DC', linestyle='solid', linewidth=0.75)
        axs[i].plot([med_ness, med_ness], [0, 0.45], c='#F5F5DC', linestyle='solid', linewidth=0.75)
        axs[i].plot([med_ess, med_ness], [0.1, 0.1], c='#F5F5DC', linestyle='solid', linewidth=0.75)

        if ms == "Organoids":
            textx, texty = ((med_ess + med_ness)/2) - 0.6, 0.8
        else:
            textx, texty = ((med_ess + med_ness)/2) - 1.1, 0.8
        arrow_dx, arrow_dy = 0, -0.14
        axs[i].text(
            textx, texty, f'$\Delta$ Median = {med_ess - med_ness:.2f}',
            horizontalalignment='center', fontsize=ANNOT_SIZE, transform=axs[i].get_xaxis_transform(), bbox=dict(pad=0.4, facecolor='white', alpha=0.5, boxstyle='Round')
            # bbox=dict(boxstyle='round', facecolor='white', alpha=1)
           )
        axs[i].arrow(textx + arrow_dx, texty + arrow_dy, ((med_ess + med_ness) / 2) - textx, (0.02 - arrow_dy) - (texty + arrow_dy), head_length=0,
                     linewidth=0.5,
                     transform=axs[i].get_xaxis_transform())

    axs[-1].set_xlabel('Mean Naive Gene Effect')

    plt.savefig(os.path.join(FIGURE_DIR, 'Fig_1g_average_nnmd_components_4cat.pdf'))
    pd.concat(ngs_controls_longform, axis=0).reset_index().rename({
        'level_0': 'Group'
    }, axis=1).drop('level_1', axis=1).to_csv(os.path.join(PROCESSED_DIR, 'average_nnmd_components_4cat.csv'), index=False)
    
# figure ED4ab
    
def plot_qc_metrics(screen_qc_table):
    """
    Plot CRISPR screen quality metrics in different groups of screens. NNMD and ROC-AUC per screen are plotted.

    Args:
        screen_qc_table (pandas.DataFrame): A table of QC metrics per screen
    """
    group_rename = {
        'Cas9 Traditional': 'Cas9 / traditional',
        'Humagne Traditional': 'enAsCas12a / traditional',
        'Organoids': 'enAsCas12a / 3d organoid',
        'Spheroids': 'enAsCas12a / NextGen CNS',
    }
    qc_categories = list(group_rename.keys())
    
    # NNMDs
    plt.figure(figsize=(80 * mm, 59 * mm))
    sns.boxplot(
        screen_qc_table,
        x='Group',
        y='ScreenNNMD',
        hue='Group', palette=depmap_set_palette, dodge=False, showfliers=False, order=qc_categories, legend=False,
        boxprops={"alpha": 0.6, 'edgecolor': 'black'},
        medianprops={'linewidth': 2}
    )

    sns.stripplot(
        screen_qc_table,
        x='Group',
        y='ScreenNNMD',
        hue='Group', palette=depmap_set_palette, dodge=False, edgecolor='black', linewidth=0.5, order=qc_categories, legend=False
    )
    
    plt.xlabel('')
    plt.xticks([0, 1, 2, 3], labels=[group_rename[x.get_text()] for x in plt.xticks()[1]], rotation=90, fontdict={'size': TICK_SIZE})

    plt.ylabel('Screen NNMD', fontdict={'size': LABEL_SIZE})
    plt.gca().tick_params(axis='x', which='both',length=0, pad=3)
    plt.yticks(plt.yticks()[0][2:-2:2], labels=[y.get_text() for y in plt.yticks()[1][2:-2:2]], fontdict={'size': TICK_SIZE})
    plt.savefig(os.path.join(FIGURE_DIR, 'Fig_ED4ab_crispr_screen_nnmds.pdf'), bbox_inches='tight')
    
    # ROC-AUCs
    plt.figure(figsize=(80 * mm, 59 * mm))
    sns.boxplot(
        screen_qc_table,
        x='Group',
        y='ScreenROCAUC',
        hue='Group', palette=depmap_set_palette, dodge=False, showfliers=False, order=qc_categories, legend=False,
        boxprops={"alpha": 0.6, 'edgecolor': 'black'},
        medianprops={'linewidth': 2}
    )

    sns.stripplot(
        screen_qc_table,
        x='Group',
        y='ScreenROCAUC',
        hue='Group', palette=depmap_set_palette, dodge=False, edgecolor='black', linewidth=0.5, order=qc_categories, legend=False
    )
    
    plt.xlabel('')
    plt.xticks([0, 1, 2, 3], labels=[group_rename[x.get_text()] for x in plt.xticks()[1]], rotation=90, fontdict={'size': TICK_SIZE})

    plt.ylabel('Screen ROC-AUC', fontdict={'size': LABEL_SIZE})
    plt.gca().tick_params(axis='x', which='both',length=0, pad=3)
    plt.yticks(plt.yticks()[0][1:-1:2], labels=[y.get_text() for y in plt.yticks()[1][1:-1:2]], fontdict={'size': TICK_SIZE})
    plt.savefig(os.path.join(FIGURE_DIR, 'Fig_ED4ab_crispr_screen_rocaucs.pdf'), bbox_inches='tight')
    
# figure 1i

def plot_library_correction_examples(screen_metadata, ngs, screen_gene_effect):
    """
    Plots gene effect scores at successive stages of normalization.

    Args:
        screen_metadata (pandas.DataFrame): A table of screens to use
        ngs (pandas.DataFrame): A matrix (screen x gene) of naive gene effects
        screen_gene_effect (pandas.DataFrame): A matrix (screen x gene) of Chronos-estimated gene effects
    """
    hu_cell_line_screens = screen_metadata[
        (screen_metadata['Library'] == 'Humagne-CD') & 
        screen_metadata['IsTraditional2D']
    ].index.tolist()
    avana_cell_line_screens = screen_metadata[
        (screen_metadata['Library'] == 'Avana') & 
        screen_metadata['IsTraditional2D']
    ].index.tolist()
    
    # naive gene effects
    ngs_hu_vs_av = pd.concat([
        ngs.loc[ngs.index.str.contains('AV'), :].mean().rename('Avana Mean'),
        ngs.loc[hu_cell_line_screens, :].mean().rename('Humagne Mean')
    ], axis=1)

    # library-corrected gene effects
    aligned_hu_vs_av = pd.concat([
        screen_gene_effect.loc[screen_gene_effect.index.str.contains('AV'), :].mean().rename('Avana Mean'),
        screen_gene_effect.loc[hu_cell_line_screens, :].mean().rename('Humagne Mean')
    ], axis=1)

    # find common axis limits
    minlim = min(ngs_hu_vs_av.min().min(), aligned_hu_vs_av.min().min())
    maxlim = min(ngs_hu_vs_av.max().max(), aligned_hu_vs_av.max().max())
    limrange = maxlim - minlim
    axlims = minlim - (0.05 * limrange), maxlim + (0.05 * limrange)
    
    # naive gene effects
    plt.figure(figsize=(24 * mm, 24 * mm))
    sns.scatterplot(
        ngs_hu_vs_av,
        x='Humagne Mean',
        y='Avana Mean',
        alpha=0.3,
        s=2,
        edgecolor='white', linewidth=0.1
    )

    plt.xlim(axlims)
    plt.ylim(axlims)
    plt.title('Naive', fontdict={'size': TITLE_SIZE})
    plt.xlabel('Humagne Mean', fontdict={'size': LABEL_SIZE})
    plt.ylabel('Avana Mean', fontdict={'size': LABEL_SIZE})
    plt.axhline(0, linestyle='dotted', color='black', zorder=-1, linewidth=0.5)
    plt.axvline(0, linestyle='dotted', color='black', zorder=-1, linewidth=0.5)
    plt.axline((0, 0), slope=1, linestyle='dashed', color='black', zorder=-1, linewidth=0.5)
    plt.xticks([-5, 0], [-5, 0], fontdict={'size': TICK_SIZE})
    plt.yticks([-5, 0], [-5, 0], fontdict={'size': TICK_SIZE})
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURE_DIR, 'Fig_1i_normalization_example_naive.pdf'))

    # library-corrected gene effects
    plt.figure(figsize=(24 * mm, 24 * mm))
    sns.scatterplot(
        aligned_hu_vs_av,
        x='Humagne Mean',
        y='Avana Mean',
        alpha=0.3,
        s=2,
        edgecolor='white', linewidth=0.1
    )

    plt.xlim(axlims)
    plt.ylim(axlims)
    plt.title('Chronos', fontdict={'size': TITLE_SIZE})
    plt.xlabel('Humagne Mean', fontdict={'size': LABEL_SIZE})
    plt.ylabel('Avana Mean', fontdict={'size': LABEL_SIZE})
    plt.axhline(0, linestyle='dotted', color='black', zorder=-1, linewidth=0.5)
    plt.axvline(0, linestyle='dotted', color='black', zorder=-1, linewidth=0.5)
    plt.axline((0, 0), slope=1, linestyle='dashed', color='black', zorder=-1, linewidth=0.5)
    plt.xticks([-5, 0], [-5, 0], fontdict={'size': TICK_SIZE})
    plt.yticks([-5, 0], [-5, 0], fontdict={'size': TICK_SIZE})
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURE_DIR, 'Fig_1i_normalization_example_final.pdf'))
    
    
def main():
    """
    Generate all Project Achilles (genome-wide CRISPR screen) quality control plots.
    """
    ngs, ness, ess, screen_dict, screen_qc_table, screen_metadata, screen_gene_effect = load_data()
    
    # figure 1g
    plot_control_separations(ngs, screen_dict, ness, ess)
    
    # figure ED3ab
    plot_qc_metrics(screen_qc_table)
    
    # figure 1i
    plot_library_correction_examples(screen_metadata, ngs, screen_gene_effect)


if __name__ == "__main__":
    main()