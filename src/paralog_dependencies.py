import os
import pandas as pd
import numpy as np
from figure_utils import *
from constants import *
from data_utils import *

import seaborn as sns
import matplotlib.pyplot as plt

mpl.style.use(os.path.join(ASSETS_DIR, 'stylesheet.mplstyle'))

# figure 2d

def plot_paralog_dependency_volcano(highlight_point=None, highlight_color='#f87060', fdr_threshold=0.005, correlation_threshold=0.3):
    """
    Plot correlation resluts for paralog dependencies

    Args:
        highlight_point (str, optional): Label indicating a point highlight, formatted as "<dependency/biomarker>". Defaults to None.
        highlight_color (str, optional): String color value in hex. Defaults to '#f87060'.
        fdr_threshold (float, optional): Minimum significance value. Defaults to 0.005.
        correlation_threshold (float, optional): Minimum effect size. Defaults to 0.3.
    """
    next_gen_paralog_dependency_volcano_df = load_file('paralog_dependency_table.csv', local_dir=PROCESSED_DIR)
    next_gen_paralog_dependency_volcano_df['p_transform'] = -np.log10(next_gen_paralog_dependency_volcano_df['NextGenPValue'])
    next_gen_paralog_dependency_volcano_df['label'] = next_gen_paralog_dependency_volcano_df['Dependency'].str.split(' ').str[0] + '/' + next_gen_paralog_dependency_volcano_df['Biomarker'].str.split(' ').str[0]

    xlabel='NextGenPearsonR'
    ylabel='p_transform'
    
    # make plot
    plt.figure(figsize=(50*mm, 50*mm))
    plt.subplots_adjust(bottom=0.15)

    # plot other points first
    sns.scatterplot(next_gen_paralog_dependency_volcano_df[next_gen_paralog_dependency_volcano_df['label'] != highlight_point],
                    x=xlabel, y=ylabel,
                    hue='Significant', palette={True: sns.color_palette('Dark2')[0], False: 'tab:gray'}, hue_order=[True, False],
                    s=4,
                    alpha = 0.6, edgecolor='black', legend=False)

    # plot the highlighted point
    sns.scatterplot(next_gen_paralog_dependency_volcano_df[next_gen_paralog_dependency_volcano_df['label'] == highlight_point],
                    x=xlabel, y=ylabel,
                    s=4,
                    color=highlight_color,
                    alpha = 0.6, edgecolor='black', legend=False)

    plt.ylabel('-log10(p-value)', fontdict={'size': LABEL_SIZE})
    plt.xlabel('Pearson Correlation\nDependency vs. Paralog Expression', fontdict={'size': LABEL_SIZE})
    plt.title('Paralog Dependencies', fontdict={'size': TITLE_SIZE})
    plt.ylim(plt.ylim()[0], plt.ylim()[1] + 0.5)

    # add a significance line
    p_fdr_approx = approximate_fdr_in_p_units(next_gen_paralog_dependency_volcano_df, fdr_threshold=fdr_threshold, pcolumn='NextGenPValue', fdr_column='FDR')
    plt.axhline(-np.log10(p_fdr_approx), linestyle='dashed', color='tab:red', alpha=0.5, linewidth=0.5)
    plt.axvline(correlation_threshold, linestyle='dashed', color='tab:red', alpha=0.5, linewidth=0.5)

    # manually add some annotations
    annots = pd.concat([
        next_gen_paralog_dependency_volcano_df.query('FDR < 0.05').sort_values('FDR')
    ], axis=0).set_index('label')
    manually_annotate(plt.gca(), 'GSPT1/GSPT2', annots, 
                      xlabel, ylabel, (-0.05, .6), ha='right', va='center',
                      arrowprops=dict(arrowstyle='-', color='black', alpha=1, linewidth=0.25, shrinkA=0, shrinkB=1),
                      bbox_kws={'facecolor': 'white', 'edgecolor': 'white', 'boxstyle':'round,pad=0.05'})
    manually_annotate(plt.gca(), 'EIF1AX/EIF1AY', annots, 
                      xlabel, ylabel, (-0.05, -1), ha='right', va='center',
                      arrowprops=dict(arrowstyle='-', color='black', alpha=1, linewidth=0.25, shrinkA=0, shrinkB=1),
                      bbox_kws={'facecolor': 'white', 'edgecolor': 'white', 'boxstyle':'round,pad=0.05'})
    manually_annotate(plt.gca(), 'DNAJC19/DNAJC15', annots, 
                      xlabel, ylabel, (-0.05, -.5), ha='right', va='center',
                      arrowprops=dict(arrowstyle='-', color='black', alpha=1, linewidth=0.25, shrinkA=0, shrinkB=1),
                      bbox_kws={'facecolor': 'white', 'edgecolor': 'white', 'boxstyle':'round,pad=0.05'})
    manually_annotate(plt.gca(), 'CDS2/CDS1', annots, 
                      xlabel, ylabel, (-0.05, -0.8), ha='right', va='center',
                      arrowprops=dict(arrowstyle='-', color='black', alpha=1, linewidth=0.25, shrinkA=0, shrinkB=1),
                      bbox_kws={'facecolor': 'white', 'edgecolor': 'white', 'boxstyle':'round,pad=0.05'})
    manually_annotate(plt.gca(), 'FAM50A/FAM50B', annots, 
                      xlabel, ylabel, (-0.07, -1.3), ha='right', va='center',
                      arrowprops=dict(arrowstyle='-', color='black', alpha=1, linewidth=0.25, shrinkA=0, shrinkB=1),
                      bbox_kws={'facecolor': 'white', 'edgecolor': 'white', 'boxstyle':'round,pad=0.05'})
    manually_annotate(plt.gca(), 'ENO1/ENO2', annots, 
                      xlabel, ylabel, (-0.06, 0), ha='right', va='center',
                      arrowprops=dict(arrowstyle='-', color='black', alpha=1, linewidth=0.25, shrinkA=0, shrinkB=1),
                      color=highlight_color,
                      bbox_kws={'facecolor': 'white', 'edgecolor': 'white', 'boxstyle':'round,pad=0.05'})

    # add a legend
    handles = []
    handles.append(Line2D([], [], color=sns.color_palette('Dark2')[0], marker = 'o', linestyle='None', 
                          markeredgewidth=0.25, markeredgecolor='black', markersize=4, label = 'Significant', alpha=0.6))
    handles.append(Line2D([], [], color='tab:gray', marker = 'o', linestyle='None', 
                          markeredgewidth=0.25, markeredgecolor='black', markersize=4, label = 'Insignificant', alpha=0.6))
    plt.legend(handles=handles, handletextpad=0, prop={'size': ANNOT_SIZE}, loc='upper left')
    
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURE_DIR, f'Fig_2d_next_gen_paralog_expression_dependencies.pdf'))

# figure 2e

def plot_paralog_example(dep, feat):
    """
    Plot scatterplot of dependency vs. expression for example paralog dependencies

    Args:
        dep (str): The dependency gene name
        feat (str): The biomarker gene name
    """
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
    expr = load_full_matrix('expression')
    screen_expr = expand_model_matrix_to_screens(expr, screen_metadata.reset_index().set_index('ModelID')['ScreenID'])
    
    next_gen_paralog_dependency_volcano_df = load_file('paralog_dependency_table.csv', local_dir=PROCESSED_DIR)
    next_gen_paralog_dependency_volcano_df['label'] = next_gen_paralog_dependency_volcano_df['Dependency'].str.split(' ').str[0] + '/' + next_gen_paralog_dependency_volcano_df['Biomarker'].str.split(' ').str[0]
    
    # map the gene names to identifiers
    feat_id = search_gene(feat)
    dep_id = search_gene(dep)
    
    # create tables for plotting example paralog relationships
    def construct_example_df(dependency, biomarker, screen_metadata_subset):
        """
        Reformat dependency and biomarker tables to facilitate making the scatterplot.

        Args:
            dependency (str): The name of the dependency gene
            biomarker (str): The name of the biomarker gene (identical to `dependency` for expression addictions)
            screen_metadata_subset (pandas.DataFrame): A table of screens to include values for

        Returns:
            pandas.DataFrame: A table of paired dependency-expression values per sample
        """
        example_df = pd.concat([
            pd.Series(dependency, index=screen_metadata_subset.index.tolist(), name='Dependency'),
            screen_gene_effect.reindex(index=screen_metadata_subset.index.tolist()).loc[:, dependency].rename('DependencyGeneEffect'),
            pd.Series(biomarker, index=screen_metadata_subset.index.tolist(), name='Biomarker'),
            screen_expr.reindex(index=screen_metadata_subset.index.tolist()).loc[:, biomarker].rename('BiomarkerExpression'),
            screen_metadata_subset['OncotreeLineage']
        ], axis=1)
        return example_df 
    next_gen_example_df = construct_example_df(dep_id, feat_id, next_gen_screen_metadata)
    adherent_example_df = construct_example_df(dep_id, feat_id, adherent_lin_matched_screen_metadata)
    
    # make pairs of regression plots on top of scatter plots
    fig, axs = plt.subplots(2, 1, figsize=(30 * mm, 50 * mm), sharey=True, gridspec_kw={'hspace': 0.1})
    plt.subplots_adjust(left=0.25, top=0.89, bottom=0.2)
    
    sns.regplot(
        next_gen_example_df,
        x='BiomarkerExpression',
        y='DependencyGeneEffect',
        scatter=False, color='gray', line_kws={'linewidth': 1},
        ax=axs[0]
    )
    sns.scatterplot(
        next_gen_example_df,
        x='BiomarkerExpression',
        y='DependencyGeneEffect',
        hue='OncotreeLineage', palette=lineage_cmap, legend=False,
        marker='s', s=4,
        ax=axs[0]
    )
    
    sns.regplot(
        adherent_example_df,
        x='BiomarkerExpression',
        y='DependencyGeneEffect',
        scatter=False, color='gray', line_kws={'linewidth': 1},
        ax=axs[1]
    )
    sns.scatterplot(
        adherent_example_df,
        x='BiomarkerExpression',
        y='DependencyGeneEffect',
        hue='OncotreeLineage', palette=lineage_cmap, legend=False,
        marker='o', s=4,
        ax=axs[1]
    )

    # add axes
    axs[0].axhline(0, color='black', linestyle='dotted', linewidth=0.5)
    axs[0].axvline(0, color='black', linestyle='dotted', linewidth=0.5)
    axs[1].axhline(0, color='black', linestyle='dotted', linewidth=0.5)
    axs[1].axvline(0, color='black', linestyle='dotted', linewidth=0.5)

    # label axes
    axs[0].set_ylabel('')
    axs[1].set_ylabel('')
    fig.supylabel(f'          {dep} Dependency', fontsize=LABEL_SIZE, x=0.03)
    axs[0].set_yticks([-2.0, -1.5, -1, -0.5, 0, 0.5], [-2.0, -1.5, -1.0, -0.5, 0.0, 0.5], fontdict={'size': TICK_SIZE})
    axs[1].set_yticks([-2.0, -1.5, -1, -0.5, 0, 0.5], [-2.0, -1.5, -1.0, -0.5, 0.0, 0.5], fontdict={'size': TICK_SIZE})

    axs[0].set_xlabel('')
    axs[1].set_xlabel(f'{feat} Expression', fontdict={'size': LABEL_SIZE})
    axs[0].set_xticks([0, 2, 4, 6, 8, 10], [], fontdict={'size': TICK_SIZE})
    axs[1].set_xticks([0, 2, 4, 6, 8, 10], [0, 2, 4, 6, 8, 10], fontdict={'size': TICK_SIZE})

    # set common axis limits
    minx = np.min([axs[0].get_xlim()[0], axs[1].get_xlim()[0]])
    maxx = np.max([axs[0].get_xlim()[1], axs[1].get_xlim()[1]])
    axs[0].set_xlim(minx, maxx)
    axs[1].set_xlim(minx, maxx)
    axs[0].set_ylim(-2.4, 0.7)
    axs[1].set_ylim(-2.4, 0.7)

    # label each figure with the dataset used
    axs[0].text(0.95, 0.88, 'NextGen', fontdict={'fontsize': ANNOT_SIZE, 'color': 'black'}, horizontalalignment='right', transform=axs[0].transAxes)
    axs[0].text(0.97, 0.05, 'r = {:.2f}'.format(next_gen_paralog_dependency_volcano_df.set_index('label').loc[dep + '/' + feat, 'NextGenPearsonR']), 
                fontdict={'fontsize': ANNOT_SIZE}, ha='right', transform=axs[0].transAxes)
    axs[1].text(0.95, 0.88, 'Traditional', fontdict={'fontsize': ANNOT_SIZE, 'color': 'black'}, horizontalalignment='right', transform=axs[1].transAxes)
    axs[1].text(0.97, 0.05, 'r = {:.2f}'.format(next_gen_paralog_dependency_volcano_df.set_index('label').loc[dep + '/' + feat, 'TraditionalPearsonR']), 
                fontdict={'fontsize': ANNOT_SIZE}, ha='right', transform=axs[1].transAxes)

    plt.savefig(os.path.join(FIGURE_DIR, f'Fig_2e_{dep}_{feat}_paralog_example.pdf'))
    
    example_summary_df = pd.concat([
        next_gen_example_df.assign(ScreenSet='NextGen'),
        adherent_example_df.assign(ScreenSet='Traditional2D')
    ], axis=0).rename_axis('ScreenID')
    example_summary_df.to_csv(os.path.join(PROCESSED_DIR, 'paralog_example_table.csv'))


def main():
    """
    Generate plots related to dependencies associated with expression of their paralog
    """
    # figure 2d
    plot_paralog_dependency_volcano(highlight_point="ENO1/ENO2", highlight_color='#f87060', fdr_threshold=0.005, correlation_threshold=0.3)
    
    # figure 2e
    plot_paralog_example('ENO1', 'ENO2')
    

if __name__ == "__main__":
    main()