import os
import pandas as pd
import numpy as np
from figure_utils import *
from constants import *
from data_utils import *

import seaborn as sns
import matplotlib.pyplot as plt

mpl.style.use(os.path.join(ASSETS_DIR, 'stylesheet.mplstyle'))

# figure 2f

def plot_onc_tsg_volcano(fdr_threshold=0.05, mean_difference_threshold=0.3):
    next_gen_onc_tsg_dependency_volcano_df = load_file('onc_tsg_dependency_table.csv', local_dir=PROCESSED_DIR)
    
    next_gen_onc_tsg_dependency_volcano_df['-log10(MWU p)'] = -np.log10(next_gen_onc_tsg_dependency_volcano_df['NextGenMannWhitneyPValue'])
    next_gen_onc_tsg_dependency_volcano_df['label'] = next_gen_onc_tsg_dependency_volcano_df['Dependency'].str.split(' ').str[0] + '/' + next_gen_onc_tsg_dependency_volcano_df['Biomarker'].str.split(' ').str[0]
    
    xlabel='NextGenMeanDifference'
    ylabel='-log10(MWU p)'
    
    plt.figure(figsize=(50 * mm, 50 * mm))
    plt.subplots_adjust(left=0.18, bottom=0.22, right=0.95, top=0.84)

    pairs_to_highlight = ["KRAS/KRAS", "SCD/KRAS"]

    sns.scatterplot(
        pd.concat([
            next_gen_onc_tsg_dependency_volcano_df[
                (next_gen_onc_tsg_dependency_volcano_df['HitClass'] != 'Insignificant') |
                (next_gen_onc_tsg_dependency_volcano_df[xlabel].abs() >= mean_difference_threshold)
            ],
            next_gen_onc_tsg_dependency_volcano_df[
                (next_gen_onc_tsg_dependency_volcano_df['HitClass'] == 'Insignificant') & 
                (next_gen_onc_tsg_dependency_volcano_df[xlabel].abs() < mean_difference_threshold)
            ].sample(frac=0.4) # reduce the density of insignificant points for plotting purposes
        ], axis=0).sample(frac=1).sort_values('HitClass', ascending=False),
        x=xlabel, y=ylabel,
        hue='HitClass', palette=onc_gof_tsg_lof_cmap,
        s=8,
        alpha=0.6, edgecolor='white',
    )
    
    p_fdr_approx = approximate_fdr_in_p_units(next_gen_onc_tsg_dependency_volcano_df, fdr_threshold=fdr_threshold, pcolumn='NextGenMannWhitneyPValue', fdr_column='FDR')
    plt.axhline(-np.log10(p_fdr_approx), linestyle='dashed', color='tab:red', alpha=0.5, linewidth=0.5)
    plt.axvline(-mean_difference_threshold, linestyle='dashed', color='tab:red', alpha=0.5, linewidth=0.5)

    annots = next_gen_onc_tsg_dependency_volcano_df[next_gen_onc_tsg_dependency_volcano_df['label'].isin(pairs_to_highlight)].set_index('label')

    manually_annotate(plt.gca(), "KRAS/KRAS", annots, 
                      xlabel, ylabel, (0.6, 1.8), ha='center', va='center',
                      arrowprops=dict(arrowstyle='-', color='black', alpha=1, linewidth=0.25, shrinkA=0, shrinkB=1),
                      bbox_kws={'facecolor': 'white', 'edgecolor': 'white', 'boxstyle':'round,pad=0.1', 'alpha':0.5, 'linewidth':0})
    manually_annotate(plt.gca(), "SCD/KRAS", annots, 
                      xlabel, ylabel, (-0.4, 1), ha='right', va='center',
                      arrowprops=dict(arrowstyle='-', color='black', alpha=1, linewidth=0.25, shrinkA=0, shrinkB=1),
                      bbox_kws={'facecolor': 'white', 'edgecolor': 'white', 'boxstyle':'round,pad=0.1', 'alpha':0.5, 'linewidth':0})

    plt.ylabel('-log10(p-value)', fontdict={'size': LABEL_SIZE})
    plt.xlabel('Mean Dependency Difference\n(Altered - Unaltered)', fontdict={'size': LABEL_SIZE})
    plt.title('Oncogene GOF and\nTumor Suppressor LOF', fontdict={'size': TITLE_SIZE})
    plt.legend(loc='upper right', handletextpad=0, prop={'size': ANNOT_SIZE})

    plt.savefig(os.path.join(FIGURE_DIR, 'Fig_2f_next_gen_onc_tsg_dependency_volcano_plot.pdf'))
    
def create_paired_violin_plot(label, onc_gof_altered_mtx, tsg_lof_altered_mtx, screen_gene_effect, screen_metadata, next_gen_screen_list, adherent_screen_list, next_gen_onc_tsg_test_df, ax_map_panel):
    next_gen_entry = next_gen_onc_tsg_test_df.set_index('label').loc[label]
    alteration_type = next_gen_entry.loc['BiomarkerAlterationType']
    dependency = next_gen_entry.loc['Dependency']
    biomarker = next_gen_entry.loc['Biomarker']
    next_gen_p = next_gen_entry.loc['NextGenMannWhitneyPValue']
    adherent_p = next_gen_entry.loc['TraditionalMannWhitneyPValue']
    
    if alteration_type == "Oncogene":
        groups = pd.DataFrame(onc_gof_altered_mtx[biomarker]).join(screen_metadata['OncotreeLineage'])
    else:
        groups = pd.DataFrame(tsg_lof_altered_mtx[biomarker]).join(screen_metadata['OncotreeLineage'])
        
    if "APC" in biomarker:
        next_gen_order = (
            screen_metadata[
                screen_metadata.index.isin(next_gen_screen_list) & 
                (screen_metadata['OncotreeLineage'] != "Bowel")
            ].index.tolist() +
            screen_metadata[
                screen_metadata.index.isin(next_gen_screen_list) & 
                (screen_metadata['OncotreeLineage'] == "Bowel")
            ].index.tolist()
        )
        adherent_order = (
            screen_metadata[
                screen_metadata.index.isin(adherent_screen_list) & 
                (screen_metadata['OncotreeLineage'] != "Bowel")
            ].index.tolist() +
            screen_metadata[
                screen_metadata.index.isin(adherent_screen_list) & 
                (screen_metadata['OncotreeLineage'] == "Bowel")
            ].index.tolist()
        )
    else:
        next_gen_order = next_gen_screen_list
        adherent_order = adherent_screen_list
        
    # plot NextGen
    make_gene_effect_distribution_plot(
        dependency,
        x=biomarker,
        box_hue=biomarker, box_palette={True: 'tab:gray', False: 'tab:gray'},
        box_alpha=0.5,
        point_hue='OncotreeLineage', point_palette=lineage_cmap,
        mtx= screen_gene_effect,
        sample_subset = next_gen_order,
        grouping_annotation = groups,
        plot_type='boxviolin',
        dodge=False,
        point_kwargs={'s': 2, 'legend': False, 'marker':'s', 'edgecolor': 'black', 'linewidth':0.25, 'jitter': 0.15},
        box_kwargs={'inner': None, 'legend': False},
        ax=ax_map_panel[(0, 0)]
    )
    
    # plot adherent
    make_gene_effect_distribution_plot(
        dependency,
        x=biomarker,
        box_hue=biomarker, box_palette={True: 'tab:gray', False: 'tab:gray'},
        box_alpha=0.5,
        point_hue='OncotreeLineage', point_palette=lineage_cmap,
        mtx= screen_gene_effect,
        sample_subset = adherent_order,
        grouping_annotation = groups,
        plot_type='boxviolin',
        dodge=False,
        point_kwargs={'s': 2, 'legend': False, 'marker':'o', 'edgecolor': 'black', 'linewidth':0.25, 'jitter': 0.15},
        box_kwargs={'inner': None, 'legend': False},
        ax=ax_map_panel[(0, 1)]
    )
    
    # Align the axes
    
    nymin, nymax = ax_map_panel[(0, 0)].get_ylim()
    aymin, aymax = ax_map_panel[(0, 1)].get_ylim()
    ymin = min(nymin, aymin)
    ymax = max(nymax, aymax)
    ybar = ymin + (0.02 * (ymax - ymin))
    ybar_text = ymin + (0.04 * (ymax - ymin))
    ybar_cap = ymin + (0.06 * (ymax - ymin))
    extend_ymin = ymin - (0.1 * (ymax - ymin))

    ax_map_panel[(0, 0)].set_ylim(extend_ymin, ymax)
    ax_map_panel[(0, 1)].set_ylim(extend_ymin, ymax)

    # Add labels
    ax_map_panel[(0, 0)].set_ylabel(dependency.split(' ')[0] + ' Dependency', fontdict={'size': LABEL_SIZE})
    ax_map_panel[(0, 1)].set_ylabel('')
    ax_map_panel[(0, 0)].set_yticks(ax_map_panel[(0, 0)].get_yticks()[1:-1], labels=ax_map_panel[(0, 0)].get_yticklabels()[1:-1], fontdict={'size': TICK_SIZE})
    ax_map_panel[(0, 1)].set_yticks(ax_map_panel[(0, 1)].get_yticks()[1:-1], labels=[], fontdict={'size': TICK_SIZE})

    ax_map_panel[(0, 0)].set_xticks([0, 1], labels=['Unaltered', 'Altered'], fontdict={'size': TICK_SIZE}, rotation=20)
    ax_map_panel[(0, 1)].set_xticks([0, 1], labels=['Unaltered', 'Altered'], fontdict={'size': TICK_SIZE}, rotation=20)
    ax_map_panel[(0, 0)].set_xlabel(biomarker.split(' ')[0] + ' Status', x=1.05, fontdict={'size': LABEL_SIZE})
    ax_map_panel[(0, 1)].set_xlabel('')

    ax_map_panel[(0, 0)].set_title('NextGen', fontdict={'size': TITLE_SIZE})
    ax_map_panel[(0, 1)].set_title('Traditional', fontdict={'size': TITLE_SIZE})

    # Add reference lines and significance
    ax_map_panel[(0, 0)].axhline(0, color='black', linestyle='dotted', linewidth=0.5)
    ax_map_panel[(0, 1)].axhline(0, color='black', linestyle='dotted', linewidth=0.5)

    ax_map_panel[(0, 0)].plot([0, 1], [ybar, ybar], color='black', linewidth=0.75)
    ax_map_panel[(0, 0)].plot([0, 0], [ybar, ybar_cap], color='black', linewidth=0.75)
    ax_map_panel[(0, 0)].plot([1, 1], [ybar, ybar_cap], color='black', linewidth=0.75)
    ax_map_panel[(0, 0)].text(
        0.5, ybar_text, format_significance(next_gen_p),
        horizontalalignment='center', verticalalignment='bottom', 
        fontsize=ANNOT_SIZE, transform=ax_map_panel[(0, 0)].transData
    )

    ax_map_panel[(0, 1)].plot([0, 1], [ybar, ybar], color='black', linewidth=0.75)
    ax_map_panel[(0, 1)].plot([0, 0], [ybar, ybar_cap], color='black', linewidth=0.75)
    ax_map_panel[(0, 1)].plot([1, 1], [ybar, ybar_cap], color='black', linewidth=0.75)
    ax_map_panel[(0, 1)].text(
        0.5, ybar_text, format_significance(adherent_p),
        horizontalalignment='center', verticalalignment='bottom', 
        fontsize=ANNOT_SIZE, transform=ax_map_panel[(0, 1)].transData
    )
    
    all_screens = next_gen_screen_list + adherent_screen_list
    example_summary_df = pd.concat([
        pd.Series(dependency, index=all_screens, name='Dependency'),
        screen_gene_effect.loc[all_screens, dependency].rename('DependencyGeneEffect'),
        pd.Series(biomarker, index=all_screens, name='Biomarker'),
        pd.Series(alteration_type, index=all_screens, name='BiomarkerAlterationType').replace({'Oncogene': 'GoF', 'Tumor Suppressor': 'LoF'}),
        groups[biomarker].rename('Altered'),
        screen_metadata.loc[all_screens, 'OncotreeLineage'],
        pd.concat([pd.Series('NextGen', index=next_gen_screen_list), pd.Series('Traditional2D', index=adherent_screen_list)]).rename('ScreenSet')
    ], axis=1).reset_index().rename({'index': 'ScreenID'}, axis=1)
    return example_summary_df.dropna()    

# figure 2g

def plot_example_panel(examples, 
                       figsize=(80 * mm, 50 * mm), 
                       column_ratios=[1, 1], column_widthspace=0.5):
    fig = plt.figure(figsize=figsize)
    plt.subplots_adjust(left=0.15, bottom=0.2, right=0.95, top=0.88)
    
    # two columns
    outer = gridspec.GridSpec(1, 2, width_ratios=column_ratios, wspace=column_widthspace)
    # left column, two rows
    gsl = gridspec.GridSpecFromSubplotSpec(subplot_spec=outer[0], nrows=1, ncols=2, wspace=0.1)
    # right column, four rows
    gsr = gridspec.GridSpecFromSubplotSpec(subplot_spec=outer[1], nrows=1, ncols=2, wspace=0.1)

    ax_map = dict()
    # left
    ax_map[0] = dict()
    ax_map[0][(0, 0)] = fig.add_subplot(gsl[0])
    ax_map[0][(0, 1)] = fig.add_subplot(gsl[1])
    # right
    ax_map[1] = dict()
    ax_map[1][(0, 0)] = fig.add_subplot(gsr[0])
    ax_map[1][(0, 1)] = fig.add_subplot(gsr[1])

    screen_gene_effect = load_file('screen_gene_effect.csv', index_col=0)
    screen_metadata = load_file('screen_metadata.csv', index_col=0)
    next_gen_onc_tsg_dependency_volcano_df = load_file('onc_tsg_dependency_table.csv', local_dir=PROCESSED_DIR)
    next_gen_onc_tsg_dependency_volcano_df['label'] = next_gen_onc_tsg_dependency_volcano_df['Dependency'].str.split(' ').str[0] + '/' + next_gen_onc_tsg_dependency_volcano_df['Biomarker'].str.split(' ').str[0]
    onc_gof_altered_mtx = load_file('oncogene_gof_alteration_matrix.csv', local_dir=PROCESSED_DIR, index_col=0)
    tsg_lof_altered_mtx = load_file('tumor_suppressor_lof_alteration_matrix.csv', local_dir=PROCESSED_DIR, index_col=0)

    next_gen_screens = screen_metadata[screen_metadata['IsNextGen'] & screen_metadata['ScreenType'].isin(['2DO', '3DO', '3DN', '2DN'])].index.tolist()
    next_gen_screens = list(set(next_gen_screens) & set(screen_gene_effect.index))
    adherent_screens = screen_metadata[screen_metadata['IsTraditional2D'] & screen_metadata['ScreenType'].isin(['2DS'])].index.tolist()
    adherent_screens = list(set(adherent_screens) & set(screen_gene_effect.index))

    example_summary_df = []
    for ax, ex in zip([ax_map[0], ax_map[1]], examples):
        example_df = create_paired_violin_plot(ex, onc_gof_altered_mtx, tsg_lof_altered_mtx, screen_gene_effect, screen_metadata, next_gen_screens, adherent_screens, next_gen_onc_tsg_dependency_volcano_df, ax)
        example_summary_df.append(example_df)

    plt.savefig(os.path.join(FIGURE_DIR, f'Fig_2g_onc_tsg_dependency_3_examples.pdf'))
    pd.concat(example_summary_df, axis=0).replace(lineage_replacement).to_csv(os.path.join(PROCESSED_DIR, 'onc_tsg_dependency_3_examples_table.csv'), index=False)

# figure 2h

def plot_onc_tsg_alteration_example(dep, feat):
    screen_metadata = load_file('screen_metadata.csv', index_col=0)
    next_gen_screen_metadata = screen_metadata[
        screen_metadata['IsNextGen'] & 
        screen_metadata['ScreenType'].str.contains('2DO|3DO|3DN|2DN')
    ]
    next_gen_lineages = next_gen_screen_metadata.value_counts('OncotreeLineage').index.tolist()
    # 2D lines in matched lineages
    adherent_lin_matched_screen_metadata = screen_metadata[
        screen_metadata['IsTraditional2D'] & 
        screen_metadata['OncotreeLineage'].isin(next_gen_lineages)
    ]
    
    screen_gene_effect = load_file('screen_gene_effect.csv', index_col=0)
    cn = load_full_matrix('copy_number')
    screen_cn = expand_model_matrix_to_screens(cn, screen_metadata.reset_index().set_index('ModelID')['ScreenID'])
    hotspot = load_full_matrix('hotspot')
    screen_hotspot = expand_model_matrix_to_screens(hotspot, screen_metadata.reset_index().set_index('ModelID')['ScreenID'])
    
    # map the gene names to identifiers
    dep_id = search_gene(dep)
    feat_id = search_gene(feat)
    
    # create tables for plotting example paralog relationships
    def construct_example_df(dependency, biomarker, screen_metadata_subset):
        example_df = pd.concat([
            pd.Series(dependency, index=screen_metadata_subset.index.tolist(), name='Dependency'),
            screen_gene_effect.reindex(index=screen_metadata_subset.index.tolist()).loc[:, dependency].rename('DependencyGeneEffect'),
            pd.Series(biomarker, index=screen_metadata_subset.index.tolist(), name='Biomarker'),
            screen_cn.reindex(index=screen_metadata_subset.index.tolist()).loc[:, biomarker].rename('BiomarkerCopyNumber'),
            (screen_hotspot > 0).reindex(index=screen_metadata_subset.index.tolist()).loc[:, biomarker].rename('BiomarkerMutation'),
            screen_metadata_subset['OncotreeLineage']
        ], axis=1)
        return example_df 
    next_gen_example_df = construct_example_df(dep_id, feat_id, next_gen_screen_metadata)
    adherent_example_df = construct_example_df(dep_id, feat_id, adherent_lin_matched_screen_metadata)

    next_gen_example_df['log(BiomarkerCopyNumber)'] = np.log2(next_gen_example_df['BiomarkerCopyNumber'] + 1)
    next_gen_example_df.loc[(next_gen_example_df['OncotreeLineage'] != 'Esophagus/Stomach') & (next_gen_example_df['BiomarkerMutation']), 'Group'] = f'Other:{feat}mut'
    next_gen_example_df.loc[(next_gen_example_df['OncotreeLineage'] != 'Esophagus/Stomach') & (~next_gen_example_df['BiomarkerMutation']), 'Group'] = f'Other:{feat}wt'
    next_gen_example_df.loc[(next_gen_example_df['OncotreeLineage'] == 'Esophagus/Stomach') & (next_gen_example_df['BiomarkerMutation']), 'Group'] = f'Eso/Stm:{feat}mut'
    next_gen_example_df.loc[(next_gen_example_df['OncotreeLineage'] == 'Esophagus/Stomach') & (~next_gen_example_df['BiomarkerMutation']), 'Group'] = f'Eso/Stm:{feat}wt'
    next_gen_r = next_gen_example_df[['DependencyGeneEffect', 'log(BiomarkerCopyNumber)']].dropna().corr().iloc[0, 1]
    
    # make pairs of regression plots on top of scatter plots
    fig = plt.figure(figsize=(50 * mm, 50 * mm))
    ax = plt.gca()
    plt.subplots_adjust(left=0.15, top=0.95, bottom=0.15, right=0.95)
    
    sns.regplot(
        next_gen_example_df,
        x='log(BiomarkerCopyNumber)',
        y='DependencyGeneEffect',
        scatter=False, color='gray', line_kws={'linewidth': 1},
        ax=ax
    )
    sns.scatterplot(
        next_gen_example_df,
        x='log(BiomarkerCopyNumber)',
        y='DependencyGeneEffect',
        hue='Group', palette={f'Other:{feat}wt': 'tab:gray', 
                              f'Other:{feat}mut': onc_gof_tsg_lof_cmap['Driven by alteration'], #'#95005c', 
                              f'Eso/Stm:{feat}wt': lineage_cmap['Esophagus/Stomach']}, legend=False,
        marker='s', s=4,
        ax=ax,
    )
    plt.ylabel(f'{dep} Dependency', fontsize=LABEL_SIZE)
    plt.xlabel(f'{feat} CN (log2[ploidy + 1])', fontsize=LABEL_SIZE)
    plt.yticks([-2, -1, 0], [-2, -1, 0], fontdict={'size': TICK_SIZE})
    plt.xticks([2, 4, 6], [2, 4, 6], fontdict={'size': TICK_SIZE})
    ax.axhline(0, color='black', linestyle='dotted', linewidth=0.5)
    ax.text(0.05, 0.05, 'r = {:.2f}'.format(next_gen_r), fontdict={'fontsize': ANNOT_SIZE}, ha='left', transform=ax.transAxes)

    handles = [
        Line2D([], [], color=lineage_cmap['Esophagus/Stomach'], marker = 's', linestyle='None', markersize=4, label = f'Eso/Stm : {feat} wt', markeredgewidth=0),
        Line2D([], [], color=onc_gof_tsg_lof_cmap['Driven by alteration'], marker = 's', linestyle='None', markersize=4, label = f'Other : {feat} mut', markeredgewidth=0),
        Line2D([], [], color='tab:gray', marker = 's', linestyle='None', markersize=4, label = f'Other : {feat} wt', markeredgewidth=0)
    ]
    plt.legend(handles=handles, handletextpad=0, prop={'size': ANNOT_SIZE}, loc='upper right')
    
    plt.savefig(os.path.join(FIGURE_DIR, f'Fig_2h_{dep}_{feat}_cn_alteration_example.pdf'))
    
    example_summary_df = pd.concat([
        next_gen_example_df.assign(ScreenSet='NextGen'),
        adherent_example_df.assign(ScreenSet='Traditional2D')
    ], axis=0).rename_axis('ScreenID')
    example_summary_df.to_csv(os.path.join(PROCESSED_DIR, f'{dep}_{feat}_cn_alteration_example_table.csv'))
    
    
def main():
    # figure 2f
    plot_onc_tsg_volcano(fdr_threshold=0.05, mean_difference_threshold=0.3)
    
    # figure 2g
    plot_example_panel(["KRAS/KRAS", "SCD/KRAS"])

    # figure 2h
    plot_onc_tsg_alteration_example('SCD', 'KRAS')
    

if __name__ == "__main__":
    main()
