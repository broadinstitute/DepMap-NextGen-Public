import os
import pandas as pd
import numpy as np
from figure_utils import *
from constants import *
from data_utils import *

import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.gridspec as gridspec
from matplotlib.colors import TwoSlopeNorm

mpl.style.use(os.path.join(ASSETS_DIR, 'stylesheet.mplstyle'))


def load_data():
    all_celligner_samples = load_file('celligner_coordinates_w_hcmi.csv', index_col=0)
    omics_models_meta = load_file('model_metadata.csv', index_col=0)
    screen_metadata = load_file('screen_metadata.csv', index_col=0)
    all_mmp_scores = load_file('celligner_geneset_scores.csv', local_dir=PROCESSED_DIR, index_col=0)
    celligner_top_nns = load_file('celligner_nearest_neighbor_table.csv', local_dir=PROCESSED_DIR)
    
    return all_celligner_samples, omics_models_meta, screen_metadata, all_mmp_scores, celligner_top_nns

def force_categories(series, category_list):
    series_remapped = series.apply(lambda x: x if x in category_list else 'Other')
    return series_remapped

# figure 1e

def global_celligner(lineages, all_celligner_samples, omics_models_meta, screen_metadata):
    # subset to the tcga table
    all_tcga_tumors_celligner_table = all_celligner_samples[
        all_celligner_samples["type"] == "TCGA+ tumor"
    ]
    next_gen_lineage_matched_tumors_celligner_table = all_tcga_tumors_celligner_table[
        all_tcga_tumors_celligner_table['lineage'].isin(lineages)
    ]

    # subset to the met500 table
    all_met500_tumors_celligner_table = all_celligner_samples[
        all_celligner_samples["type"] == "MET500 tumor"
    ]
    next_gen_lineage_matched_mets_celligner_table = all_met500_tumors_celligner_table[
        all_met500_tumors_celligner_table['lineage'].isin(lineages)
    ]

    # subset to the hcmi table
    all_hcmi_tumors_celligner_table = all_celligner_samples[
        all_celligner_samples["type"] == "HCMI tumor"
    ]
    next_gen_lineage_matched_hcmi_celligner_table = all_hcmi_tumors_celligner_table[
        all_hcmi_tumors_celligner_table['lineage'].isin(lineages)
    ]

    # combine to an all-tumor table
    all_tumors_celligner_table = pd.concat([
        all_tcga_tumors_celligner_table,
        all_met500_tumors_celligner_table#, 
        # all_hcmi_tumors_celligner_table
    ], axis=0)

    # subset to the depmap models available for use first
    all_depmap_models_celligner_table = all_celligner_samples[
        all_celligner_samples["ModelID"].isin(omics_models_meta.index.tolist()) &
        ~all_celligner_samples.index.isin(screen_metadata[screen_metadata['ScreenType'] == '2DO']['ModelID'].tolist())
    ]

    nextgen_dict = omics_models_meta["IsNextGen"].to_dict()
    traditional_dict = omics_models_meta["IsTraditional2D"].to_dict()
    adherent_dict = omics_models_meta["IsAdherent2D"].to_dict()
    
    all_depmap_models_celligner_table["IsNextGen"] = all_depmap_models_celligner_table["ModelID"].map(nextgen_dict)
    all_depmap_models_celligner_table["IsTraditional2D"] = all_depmap_models_celligner_table["ModelID"].map(traditional_dict)
    all_depmap_models_celligner_table["IsAdherent2D"] = all_depmap_models_celligner_table["ModelID"].map(adherent_dict)

    # from those models, separate them into 3d and 2d
    next_gen_depmap_models_celligner_table = all_depmap_models_celligner_table[
        all_depmap_models_celligner_table['IsNextGen']
    ]
    cell_line_depmap_models_celligner_table = all_depmap_models_celligner_table[
        all_depmap_models_celligner_table['IsTraditional2D']
    ]
    
    tcga_marker = 'P'
    tcga_pointsize = 2
    tcga_alpha = 0.3
    met500_marker = 'X'
    met500_pointsize = 2
    met500_alpha = 0.3
    hcmi_marker = 'o'
    hcmi_pointsize = 2
    hcmi_alpha = 0.3
    nextgen_marker = 's'
    nextgen_pointsize = 4
    nextgen_alpha = 1

    plt.figure(figsize=(68*mm, 36*mm))
    plt.subplots_adjust(left=0.06, bottom=0.04, top=0.98, right=1)

    # plot background TCGA
    celligner_tcga_df = all_tcga_tumors_celligner_table.copy()
    celligner_tcga_df['lineage'] = force_categories(celligner_tcga_df['lineage'], lineages)
    sns.scatterplot(celligner_tcga_df, x='umap1', y='umap2', s=tcga_pointsize, linewidth=0.05,
                    hue='lineage', palette=lineage_cmap,
                    marker=tcga_marker, alpha=tcga_alpha, zorder=0)

    # plot background Met500
    celligner_met500_df = all_met500_tumors_celligner_table.copy()
    celligner_met500_df['lineage'] = force_categories(celligner_met500_df['lineage'], lineages)
    sns.scatterplot(celligner_met500_df, x='umap1', y='umap2', s=met500_pointsize, linewidth=0.05,
                    hue='lineage', palette=lineage_cmap, 
                    marker=met500_marker, alpha=met500_alpha, zorder=0)

    # plot foreground NextGen DepMap models
    celligner_nextgen_df = next_gen_depmap_models_celligner_table.copy()
    celligner_nextgen_df['lineage'] = force_categories(celligner_nextgen_df['lineage'], lineages)
    sns.scatterplot(celligner_nextgen_df, x='umap1', y='umap2', s=nextgen_pointsize, 
                    hue='lineage', palette=lineage_cmap, edgecolor='black',
                    marker=nextgen_marker, alpha=nextgen_alpha, zorder=1)

    # add legend
    handles = []
    handles.append(Line2D([], [], color='black', marker = nextgen_marker, linestyle='None', markeredgewidth=1,
                          markersize=4, label = 'NextGen'))
    handles.append(Line2D([], [], color='black', marker = tcga_marker, linestyle='None', markeredgewidth=1,
                          markersize=4, label = 'TCGA+'))
    handles.append(Line2D([], [], color='black', marker = met500_marker, linestyle='None', markeredgewidth=1,
                          markersize=4, label = 'Met500'))
    # handles.append(Line2D([], [], color='black', marker = hcmi_marker, linestyle='None', markeredgewidth=1,
    #                       markersize=4, label = 'HCMI'))

    leg = plt.legend(handles = handles, bbox_to_anchor=(0, 1), loc='upper left', prop={'size': ANNOT_SIZE})
    plt.xlabel('')
    plt.ylabel('')
    plt.xticks([])
    plt.yticks([])
    plt.gca().spines[['top', 'bottom', 'left', 'right']].set_visible(False)

    # add axes
    xmin, xmax = plt.xlim()
    ymin, ymax = plt.ylim()
    plt.arrow(x=xmin+1.5, y=ymin+1, dx=0.1 * (xmax - xmin), dy=0, head_width=1.5, head_length=0.5, overhang=1, linewidth=0.5)
    plt.arrow(x=xmin+1.5, y=ymin+1, dx=0, dy=0.1 * (xmax - xmin), head_width=1.5, head_length=0.5, overhang=1, linewidth=0.5)
    plt.text(0.06, 0.1, 'UMAP2', fontdict={'fontsize': LABEL_SIZE}, transform=plt.gca().transAxes, horizontalalignment='right')
    plt.text(0.065, -0.01, 'UMAP1', fontdict={'fontsize': LABEL_SIZE}, transform=plt.gca().transAxes)

    plt.savefig(os.path.join(FIGURE_DIR, 'Fig_1e_celligner_next_gen_tcga_met500_hcmi.pdf'), dpi=1200)
    
# figure 4c

def pdac_celligner(all_celligner_samples, omics_models_meta, all_mmp_scores):
    
    # reformat celligner_samples to index by ModelID if it exists
    all_celligner_samples_other = all_celligner_samples[
        all_celligner_samples["ModelID"].isna()
    ]
    all_celligner_samples_depmap = all_celligner_samples[
        ~all_celligner_samples["ModelID"].isna()
    ]
    all_celligner_samples_depmap.index = all_celligner_samples_depmap.ModelID
    all_celligner_samples = pd.concat(
        [all_celligner_samples_depmap, all_celligner_samples_other]
    )

    adherent2d_models = omics_models_meta[
        omics_models_meta["IsAdherent2D"]
    ].index.tolist()
    nextgen_models = omics_models_meta[
        omics_models_meta["IsNextGen"]
    ].index.tolist()

    # filter models to use
    celligner_full_background = pd.concat([
        all_celligner_samples[
            (all_celligner_samples['type'] == "DepMap Model") &
            (all_celligner_samples["ModelID"].isin(adherent2d_models))
        ].assign(sample_set = 'DepMap 2D'),
        all_celligner_samples[
            (all_celligner_samples['type'] == "DepMap Model") &
            (all_celligner_samples["ModelID"].isin(nextgen_models))
        ].assign(sample_set = 'DepMap Organoid'),
        all_celligner_samples[
            (all_celligner_samples['type'] == "TCGA+ tumor")
        ].assign(sample_set = 'TCGA+'),
    ]).join(all_mmp_scores['PDAC-classical']).replace(lineage_replacement)

    # define foreground
    lineages_to_show = ['Pancreas', 'Esophagus/Stomach', 'Colorectal']
    celligner_full_foreground = celligner_full_background[
        celligner_full_background['lineage'].isin(lineages_to_show)
    ]

    # match the axis limits
    linxmin, linymin = celligner_full_foreground[['umap1', 'umap2']].min().values.tolist()
    linxmax, linymax = celligner_full_foreground[['umap1', 'umap2']].max().values.tolist()
    xmin = linxmin - (0.05 * (linxmax - linxmin))
    xmax = linxmax + (0.05 * (linxmax - linxmin))
    ymin = linymin - (0.05 * (linymax - linymin))
    ymax = linymax + (0.05 * (linymax - linymin))
    celligner_style_mapping = {'DepMap 2D': 'o', 'DepMap Organoid': 's', 'TCGA+': 'P', 
                               'Pancreas': 'H', 'Esophagus/Stomach': (4, 1, 45), 'Colorectal': 'd'}

    # define the coloring
    pdac_color_norm = TwoSlopeNorm(
        vcenter=np.nanquantile(celligner_full_foreground.query('sample_set.isin(["DepMap Organoid"])')['PDAC-classical'].dropna(), 0.03),
        vmin=np.nanquantile(all_mmp_scores['PDAC-classical'], 0),
        vmax=np.nanquantile(all_mmp_scores['PDAC-classical'], 0.99),
    )
    pdac_palette = sns.color_palette('viridis', as_cmap=True)
    
    # make a series of vertical plots
    fig, ax_map = prepare_gridspec(
        pd.DataFrame([[True, True]]),
        figsize = (50 * mm, 70 * mm),
        width_ratios=[24, 1], wspace=0.3,
        inner_gs_dict={
            (0, 0): {'nrows':3, 'ncols':1, 'hspace': 0.1}
        }, 
    )
    plt.subplots_adjust(top=0.95, bottom=0.1, right=0.88, left=0.1)

    for i, pair in enumerate(zip([['TCGA+'], ['DepMap Organoid'], ['DepMap 2D']], [['DepMap Organoid', 'DepMap 2D'], ['TCGA+'], ['TCGA+']])):
        foreground_subset = pair[0]
        midground_subset = pair[1]

        # plot out of lineage, opposite dataset
        background = celligner_full_background[
            ~celligner_full_background['lineage'].isin(lineages_to_show)
        ]
        sns.scatterplot(
            background.sort_values('PDAC-classical'),
            x='umap1', y='umap2', s=2,
            color='lightgray', alpha=0.4, linewidth=0.1,
            style='sample_set', markers=celligner_style_mapping,
            ax=ax_map[0][0][(i, 0)], legend=False
        )

        # plot in-lineage opposite dataset
        midground = celligner_full_foreground[
            celligner_full_foreground['lineage'].isin(lineages_to_show) &
            celligner_full_foreground['sample_set'].isin(midground_subset)
        ]
        sns.scatterplot(
            midground.sort_values('PDAC-classical'),
            x='umap1', y='umap2', s=3,
            color='tab:gray', alpha=0.4, linewidth=0.1,
            style='lineage', markers=celligner_style_mapping,
            ax=ax_map[0][0][(i, 0)], legend=False
        )

        # plot in-lineage, within dataset
        foreground = celligner_full_foreground[
            celligner_full_foreground['lineage'].isin(lineages_to_show) &
            celligner_full_foreground['sample_set'].isin(foreground_subset)
        ]
        sns.scatterplot(
            foreground.sort_values('PDAC-classical'),
            x='umap1', y='umap2', s=4,
            hue='PDAC-classical', palette=pdac_palette, hue_norm=pdac_color_norm, alpha=0.6, linewidth=0.1,
            style='lineage', markers=celligner_style_mapping, edgecolor='black',
            ax=ax_map[0][0][(i, 0)], legend=False,
        )

        # adjust the plot limits and stylings
        ax_map[0][0][(i, 0)].set_xlim(xmin, xmax)
        ax_map[0][0][(i, 0)].set_ylim(ymin, ymax)
        ax_map[0][0][(i, 0)].spines[['top', 'bottom', 'left', 'right']].set_visible(True)
        ax_map[0][0][(i, 0)].set_xticks([])
        ax_map[0][0][(i, 0)].set_yticks([])
        ax_map[0][0][(i, 0)].set_xlabel('')
        ax_map[0][0][(i, 0)].set_ylabel('')
        ax_map[0][0][(i, 0)].text(0.03, 0.05, foreground_subset[0], fontdict={'fontsize': LABEL_SIZE}, transform=ax_map[0][0][(i, 0)].transAxes)

    # add a colorbar
    sm = plt.cm.ScalarMappable(cmap=pdac_palette, norm=pdac_color_norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, cax=ax_map[0][1], orientation='vertical')
    cbar.ax.yaxis.set_label_position('left')
    cbar.set_ticks([-2, 0, 2, 4], labels=[-2, 0, 2, 4], fontsize=TICK_SIZE)
    cbar.set_label('PDAC-classical/EED Score', rotation=90, labelpad=3, fontsize=LABEL_SIZE)
    cbar.outline.set_linewidth(0.25)
    cbar.outline.set_color('black')

    # add a legend
    handles = []
    handles.append(Line2D([], [], color='black', marker = celligner_style_mapping['Pancreas'], 
                          linestyle='None', markeredgewidth=0.5, markersize=4, label = 'Pancreas'))
    handles.append(Line2D([], [], color='black', marker = celligner_style_mapping['Esophagus/Stomach'], 
                          linestyle='None', markeredgewidth=0.5, markersize=4, label = 'Esophagus/Stomach'))
    handles.append(Line2D([], [], color='black', marker = celligner_style_mapping['Colorectal'], 
                          linestyle='None', markeredgewidth=0.5, markersize=4, label = 'Colorectal'))
    ax_map[0][0][(2, 0)].legend(handles = handles, loc='lower left', bbox_to_anchor=(-0.09, -0.3), ncol=3, 
                                handletextpad=-0.3, columnspacing=0.2, prop={'size': ANNOT_SIZE})
    plt.savefig(os.path.join(FIGURE_DIR, 'Fig_4c_celligner_gavish_pdac_classical.pdf'), dpi=1200)

# figure 1f
    
def celligner_classification(next_gen_lineage_order, celligner_top_nns):
    # lineage_abbreviations = {'CNS/Brain': 'CNS/Brain', 'Esophagus/Stomach': 'Eso/Stm', 'Pancreas': 'Pancreas', 'Bowel': 'Bowel', 'Breast': 'Breast', 'Ovary/Fallopian Tube': 'Ovary', 'Prostate': 'Prostate'}
    lineage_abbreviations = {'CNS/Brain': 'CNS/Brain', 'Esophagus/Stomach': 'Eso/Stm', 'Pancreas': 'Pancreas', 'Colorectal': 'Colorectal', 'Breast': 'Breast', 'Ovary/Fallopian Tube': 'Ovary', 'Prostate': 'Prostate'}

    def generate_normalized_confusion_matrix(input_boolean_column, max_rank):
        # filter input models
        subset_table = celligner_top_nns[
            (celligner_top_nns['rank'] <= max_rank) & 
            (celligner_top_nns[input_boolean_column])
        ].replace(lineage_replacement)
        
        # pivot into matrix of counts
        confusion = subset_table.pivot_table(columns='input_lineage', index='target_lineage', values='input', aggfunc='count').fillna(0)
        total_neighbors_by_lineage = subset_table.groupby('input_lineage').count()['target']
        
        # divide by lineages to obtain fractions
        confusion_heatmap = (confusion / total_neighbors_by_lineage).reindex(
            index=next_gen_lineage_order, columns=next_gen_lineage_order
        ).fillna(0)
        
        return confusion_heatmap
        
    next_gen_confusion_heatmap = generate_normalized_confusion_matrix('inputIsNextGen', 25)
    adherent_confusion_heatmap = generate_normalized_confusion_matrix('inputIsTraditional2D', 25)
    
    # produce the paired heatmaps
    fig, ax_map = prepare_gridspec(
        pd.DataFrame([[True, True]]),
        figsize = (68 * mm, 34 * mm),
        width_ratios=[24, 1], wspace=0.2,
        inner_gs_dict={
            (0, 0): {'nrows':1, 'ncols':2, 'wspace': 0.1}
        }, 
    )
    plt.subplots_adjust(left=0.34, bottom=0.38, top=0.84)

    sns.heatmap(adherent_confusion_heatmap, vmin=0, vmax=1, cmap='Blues', ax=ax_map[0][0][(0, 1)], cbar=False,
                linecolor='gray', linewidths=0.02, xticklabels=True, yticklabels=True)

    sns.heatmap(next_gen_confusion_heatmap, vmin=0, vmax=1, cmap='Blues', ax=ax_map[0][0][(0, 0)], cbar=False,
                linecolor='gray', linewidths=0.02, xticklabels=True, yticklabels=True)
    ax_map[0][0][(0, 1)].set_ylabel('')
    ax_map[0][0][(0, 1)].set_yticks([])
    ax_map[0][0][(0, 0)].set_yticks([0.5 + i for i in range(len(next_gen_lineage_order))], next_gen_lineage_order, fontdict={'size': TICK_SIZE})
    ax_map[0][0][(0, 0)].set_xticks([0.5 + i for i in range(len(next_gen_lineage_order))], [lineage_abbreviations[lin] for lin in next_gen_lineage_order], fontdict={'size': TICK_SIZE})
    ax_map[0][0][(0, 1)].set_xticks([0.5 + i for i in range(len(next_gen_lineage_order))], [lineage_abbreviations[lin] for lin in next_gen_lineage_order], fontdict={'size': TICK_SIZE})

    ax_map[0][0][(0, 0)].set_ylabel('Tumor Lineage', fontdict={'size': LABEL_SIZE})
    ax_map[0][0][(0, 0)].set_xlabel('Model Lineage', x=1.05)
    ax_map[0][0][(0, 1)].set_xlabel('')

    ax_map[0][0][(0, 1)].set_title('Traditional', fontdict={'size': TITLE_SIZE})
    ax_map[0][0][(0, 0)].set_title('NextGen', fontdict={'size': TITLE_SIZE})

    sm = plt.cm.ScalarMappable(cmap='Blues', norm=Normalize(vmin=0, vmax=1))
    sm.set_array([])
    cbar = fig.colorbar(sm, cax=ax_map[0][1], shrink=0.1)
    cbar.ax.set_ylabel('Fraction NN Matching Lineage', rotation=90, labelpad=2, fontdict={'size': ANNOT_SIZE})
    cbar.ax.set_yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0], [0.0, 0.2, 0.4, 0.6, 0.8, 1.0], fontdict={'size': ANNOT_SIZE})
    cbar.ax.yaxis.set_label_position('left')
    cbar.outline.set_linewidth(0.25)
    plt.savefig(os.path.join(FIGURE_DIR, 'Fig_1f_celligner_classification_confusion_lineage_subset.pdf'))
    next_gen_confusion_heatmap.to_csv(os.path.join(PROCESSED_DIR, 'next_gen_celligner_classification_confusion_lineage_subset.csv'))
    adherent_confusion_heatmap.to_csv(os.path.join(PROCESSED_DIR, 'traditional_celligner_classification_confusion_lineage_subset.csv'))
    
# figure ED2b
    
def celligner_extended_classification(next_gen_lineage_order, celligner_top_nns):
    # lineage_abbreviations = {'CNS/Brain': 'CNS/Brain', 'Esophagus/Stomach': 'Eso/Stm', 'Pancreas': 'Pancreas', 'Bowel': 'Bowel', 'Breast': 'Breast', 'Ovary/Fallopian Tube': 'Ovary', 'Prostate': 'Prostate'}
    lineage_abbreviations = {'CNS/Brain': 'CNS/Brain', 'Esophagus/Stomach': 'Eso/Stm', 'Pancreas': 'Pancreas', 'Colorectal': 'Colorectal', 'Breast': 'Breast', 'Ovary/Fallopian Tube': 'Ovary', 'Prostate': 'Prostate'}

    def generate_normalized_confusion_matrix(input_boolean_column, max_rank):
        # filter input models
        subset_table = celligner_top_nns[
            (celligner_top_nns['rank'] <= max_rank) & 
            (celligner_top_nns[input_boolean_column])
        ].replace(lineage_replacement)
        
        # pivot into matrix of counts
        confusion = subset_table.pivot_table(columns='input_lineage', index='target_lineage', values='input', aggfunc='count').fillna(0)
        total_neighbors_by_lineage = subset_table.groupby('input_lineage').count()['target']
        
        # divide by lineages to obtain fractions
        confusion_heatmap = (confusion / total_neighbors_by_lineage)
        
        return confusion_heatmap
        
    adherent_full_confusion_heatmap = generate_normalized_confusion_matrix('inputIsTraditional2D', 25)
    next_gen_full_confusion_heatmap = generate_normalized_confusion_matrix('inputIsNextGen', 25)
    remaining_next_gen = sorted(list(set(next_gen_full_confusion_heatmap.columns.tolist()) - set(next_gen_lineage_order)))
    common_lineages = sorted(list(
        set(adherent_full_confusion_heatmap.index.tolist()) & set(adherent_full_confusion_heatmap.columns.tolist())
    ))
    all_next_gen_lineages = next_gen_lineage_order + remaining_next_gen
    next_gen_then_common_lineages = all_next_gen_lineages + sorted(list(set(common_lineages) - set(all_next_gen_lineages)))
    
    adherent_full_confusion_heatmap = adherent_full_confusion_heatmap.reindex(index=next_gen_then_common_lineages, columns=next_gen_lineage_order)
    next_gen_full_confusion_heatmap = next_gen_full_confusion_heatmap.reindex_like(adherent_full_confusion_heatmap).fillna(0)
    
    fig, ax_map = prepare_gridspec(
        pd.DataFrame([[True, True]]),
        figsize = (68 * mm, 75 * mm),
        width_ratios=[24, 1], wspace=0.2,
        inner_gs_dict={
            (0, 0): {'nrows':1, 'ncols':2, 'wspace': 0.1}
        }, 
    )

    plt.subplots_adjust(left=0.42, bottom=0.2, top=0.92)

    sns.heatmap(
        adherent_full_confusion_heatmap, 
        vmin=0, vmax=1, cmap='Blues', #cbar_kws={'label': 'Fraction NN Matching Lineage'},
        linecolor='gray', linewidths=0.02, xticklabels=True, yticklabels=True, 
        cbar=False,
        ax=ax_map[0][0][(0, 1)]
    )

    sns.heatmap(
        next_gen_full_confusion_heatmap, 
        vmin=0, vmax=1, cmap='Blues', #cbar_kws={'label': 'Fraction NN Matching Lineage'},
        linecolor='gray', linewidths=0.02, xticklabels=True, yticklabels=True, 
        cbar=False,
        ax=ax_map[0][0][(0, 0)]
    )

    ax_map[0][0][(0, 1)].set_title('Traditional')
    ax_map[0][0][(0, 0)].set_ylabel('Tumor Lineage')
    ax_map[0][0][(0, 0)].set_xticks([0.5 + i for i in range(len(next_gen_lineage_order))], [lineage_abbreviations[lin] for lin in next_gen_lineage_order], fontdict={'size': TICK_SIZE})

    ax_map[0][0][(0, 0)].set_title('NextGen')
    ax_map[0][0][(0, 1)].set_yticks([])
    ax_map[0][0][(0, 1)].set_xticks([0.5 + i for i in range(len(next_gen_lineage_order))], [lineage_abbreviations[lin] for lin in next_gen_lineage_order], fontdict={'size': TICK_SIZE})
    ax_map[0][0][(0, 1)].set_ylabel('')
    ax_map[0][0][(0, 0)].set_xlabel('Model Lineage', x=1.05)
    ax_map[0][0][(0, 1)].set_xlabel('')

    sm = plt.cm.ScalarMappable(cmap='Blues', norm=Normalize(vmin=0, vmax=1))
    sm.set_array([])
    cbar = fig.colorbar(sm, cax=ax_map[0][1], shrink=0.1)
    cbar.ax.set_ylabel('Fraction NN Matching Lineage', rotation=90, labelpad=2, fontdict={'size': ANNOT_SIZE})
    cbar.ax.set_yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0], [0.0, 0.2, 0.4, 0.6, 0.8, 1.0], fontdict={'size': ANNOT_SIZE})
    cbar.ax.yaxis.set_label_position('left')
    cbar.outline.set_linewidth(0.25)

    plt.savefig(os.path.join(FIGURE_DIR, 'Fig_ED2b_celligner_classification_full_confusion_comparison.pdf'))
    next_gen_full_confusion_heatmap.to_csv(os.path.join(PROCESSED_DIR, 'next_gen_celligner_classification_confusion_lineage_full.csv'))
    adherent_full_confusion_heatmap.to_csv(os.path.join(PROCESSED_DIR, 'traditional_celligner_classification_confusion_lineage_full.csv'))
    
    
def main():
    all_celligner_samples, omics_models_meta, screen_metadata, all_mmp_scores, celligner_top_nns = load_data()
    
    next_gen_screened_lineages = omics_models_meta[
        omics_models_meta['IsNextGen'] &
        (~omics_models_meta['OncotreeLineage'].isin(['Lung', 'Head and Neck', 'Other']))
    ]['OncotreeLineage'].value_counts().index.tolist()
    
    # figure 1e
    global_celligner(next_gen_screened_lineages, all_celligner_samples, omics_models_meta, screen_metadata)
    
    # figure 4c
    pdac_celligner(all_celligner_samples, omics_models_meta, all_mmp_scores)
    
    next_gen_lineage_order = [
        'CNS/Brain', 'Esophagus/Stomach', 'Pancreas', 'Colorectal', 
        'Breast', 'Ovary/Fallopian Tube', 'Prostate'
    ]
    # figure 1f
    celligner_classification(next_gen_lineage_order, celligner_top_nns)\
    
    # figure ED1c
    celligner_extended_classification(next_gen_lineage_order, celligner_top_nns)
    

if __name__ == "__main__":
    main()
