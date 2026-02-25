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

# figure 2a

def plot_expression_addiction_volcano(fdr_threshold=0.005, correlation_threshold=0.3, manual_annotate=MANUAL_ANNOTATE):
    """
    Plot correlation results for expression addiction dependencies.

    Args:
        fdr_threshold (float, optional): _description_. Defaults to 0.005.
        correlation_threshold (float, optional): _description_. Defaults to 0.3.
        manual_annotate (_type_, optional): _description_. Defaults to MANUAL_ANNOTATE.
    """
    next_gen_expr_addiction_volcano_df = load_file('expression_addiction_table.csv', local_dir=PROCESSED_DIR, index_col=0)

    next_gen_expr_addiction_volcano_df['FunctionalClass'] = pd.Categorical(next_gen_expr_addiction_volcano_df['FunctionalClass'].values.tolist(), categories=['Oncogenic TF', 'Oncogene', 'TF', 'Other'], ordered=True)
    next_gen_expr_addiction_volcano_df['p_transform'] = -np.log10(next_gen_expr_addiction_volcano_df['NextGenPValue'])
    
    xlabel='NextGenPearsonR'
    ylabel='p_transform'
    
    # make the plot
    plt.figure(figsize=(50 * mm, 50 * mm))
    plt.subplots_adjust(bottom=0.15)

    sns.scatterplot(
        next_gen_expr_addiction_volcano_df.sort_values('FunctionalClass', ascending=False),
        x=xlabel, y=ylabel,
        hue='FunctionalClass', palette=addiction_class_cmap,
        s=4,
        alpha=0.6, edgecolor='black', linewidth=0.25
    )
    plt.xlabel('Pearson Correlation\nDependency vs. Self Expression', fontdict={'size': LABEL_SIZE})
    plt.ylabel('-log10(p-value)', fontdict={'size': LABEL_SIZE})
    plt.title('Expression Addictions', fontdict={'size': TITLE_SIZE})

    # add a significance line
    p_fdr_approx = approximate_fdr_in_p_units(next_gen_expr_addiction_volcano_df, fdr_threshold=fdr_threshold, pcolumn='NextGenPValue', fdr_column='FDR')
    plt.axhline(-np.log10(p_fdr_approx), linestyle='dashed', color='tab:red', alpha=0.5, linewidth=0.5)
    plt.axvline(-correlation_threshold, linestyle='dashed', color='tab:red', alpha=0.5, linewidth=0.5)

    # annotate the top genes
    annots = next_gen_expr_addiction_volcano_df[
        (next_gen_expr_addiction_volcano_df['FDR'] < fdr_threshold) &
        ~next_gen_expr_addiction_volcano_df['FunctionalClass'].isin(['Other'])
    ].sort_values(['FDR', 'NextGenPValue']).assign(label=lambda x: x['Biomarker'].str.split(' ').str[0]).set_index('label')
    
    if manual_annotate:
        manually_annotate(plt.gca(), 'SOX2', annots, 
                          xlabel, ylabel, (0.05, 0), ha='left', va='center',
                          arrowprops=dict(arrowstyle='-', color='black', alpha=1, linewidth=0.25, shrinkA=0, shrinkB=1))
        manually_annotate(plt.gca(), 'HNF4A', annots, 
                          xlabel, ylabel, (0.4, -0.75), ha='center', va='center',
                          arrowprops=dict(arrowstyle='-', color='black', alpha=1, linewidth=0.25, shrinkA=0, shrinkB=1))
        manually_annotate(plt.gca(), 'TP63', annots, 
                          xlabel, ylabel, (0.13, -0.5), ha='left', va='center',
                          arrowprops=dict(arrowstyle='-', color='black', alpha=1, linewidth=0.25, shrinkA=0, shrinkB=1),
                          bbox_kws={'facecolor': 'white', 'edgecolor': 'white', 'boxstyle':'round,pad=0.05'})
        manually_annotate(plt.gca(), 'TCF7L2', annots, 
                          xlabel, ylabel, (0.1, -1.4), ha='left', va='center',
                          arrowprops=dict(arrowstyle='-', color='black', alpha=1, linewidth=0.25, shrinkA=0, shrinkB=1),
                          bbox_kws={'facecolor': 'white', 'edgecolor': 'white', 'boxstyle':'round,pad=0.05'})
        manually_annotate(plt.gca(), 'NFIB', annots, 
                          xlabel, ylabel, (0.3, -0.3), ha='left', va='center',
                          arrowprops=dict(arrowstyle='-', color='black', alpha=1, linewidth=0.25, shrinkA=0, shrinkB=1), 
                          bbox_kws={'facecolor': 'white', 'edgecolor': 'white', 'boxstyle':'round,pad=0.05'})
        manually_annotate(plt.gca(), 'CDK6', annots, 
                          xlabel, ylabel, (0.09, -0.7), ha='left', va='center',
                          arrowprops=dict(arrowstyle='-', color='black', alpha=1, linewidth=0.25, shrinkA=0, shrinkB=1),
                          bbox_kws={'facecolor': 'white', 'edgecolor': 'white', 'boxstyle':'round,pad=0.05'})
        manually_annotate(plt.gca(), 'CCND2', annots, 
                          xlabel, ylabel, (0.08, -1.2), ha='left', va='center',
                          arrowprops=dict(arrowstyle='-', color='black', alpha=1, linewidth=0.25, shrinkA=0, shrinkB=1),
                          bbox_kws={'facecolor': 'white', 'edgecolor': 'white', 'boxstyle':'round,pad=0.05'})
    else:
        texts = [plt.gca().text(r[xlabel], r[ylabel], i, ha='center', va='center', fontsize=ANNOT_SIZE) for i,r in annots.iterrows()]

    # add a legend
    handles = []
    for func_class in ['Oncogenic TF', 'Oncogene', 'TF', 'Other']:
        handles.append(Line2D([], [], color=addiction_class_cmap[func_class], marker = 'o', linestyle='None', markeredgewidth=0.25, markeredgecolor='black', markersize=4, label = func_class, alpha = 0.6))

    plt.legend(handles=handles, handletextpad=0, prop={'size': ANNOT_SIZE}, loc='upper right')
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURE_DIR, f'Fig_2a_next_gen_expression_addiction_volcano.pdf'))
    
# figure 2b
    
def plot_expression_addiction_enrichment_barplots():
    """
    Plot stacked barplots for enrichment of oncogenes and transcription factors within the class of expression addiction dependencies
    """
    geneset_table = load_file('geneset_table_updated.csv')
    gene_functional_class_df = geneset_table[
        geneset_table['Geneset'].isin(['Oncogenes', 'Transcription factors'])
    ].assign(indicator=True).pivot(index='cds_gene_id', columns='Geneset', values='indicator').fillna(False)
    
    next_gen_expr_addiction_volcano_df = load_file('expression_addiction_table.csv', local_dir=PROCESSED_DIR, index_col=0)
    
    addiction_class_enrichment_df = next_gen_expr_addiction_volcano_df.join(gene_functional_class_df[['Oncogenes', 'Transcription factors']]).fillna(False)
    
    # cross-tabulate significant hits against oncogene and TF annotations
    oncogene_crosstab_df = pd.crosstab(
        addiction_class_enrichment_df['Oncogenes'].replace({True: 'True', False: 'False'}), 
        addiction_class_enrichment_df['Significant'].replace({True: 'True', False: 'False'})
    )
    oncogene_crosstab_pct_df = oncogene_crosstab_df / oncogene_crosstab_df.sum()
    tf_crosstab_df = pd.crosstab(
        addiction_class_enrichment_df['Transcription factors'].replace({True: 'True', False: 'False'}), 
        addiction_class_enrichment_df['Significant'].replace({True: 'True', False: 'False'})
    )
    tf_crosstab_pct_df = tf_crosstab_df / tf_crosstab_df.sum()
    
    # make the plots
    fig, axs = plt.subplots(1, 2, figsize=(45 * mm, 50 * mm), sharey=True)
    plt.subplots_adjust(left=0.22, bottom=0.2, right=0.95, top=0.9)

    # plot the oncogenes in the left subplot, with percentages on the axes and numbers in the annotations
    LABEL_VERTICAL_ADJUST = 0.006
    axs[0].bar(['Expr. addiction', 'Other genes'], oncogene_crosstab_pct_df.loc['False', ['True', 'False']])
    axs[0].bar(['Expr. addiction', 'Other genes'], oncogene_crosstab_pct_df.loc['True', ['True', 'False']], 
               bottom=oncogene_crosstab_pct_df.loc['False', ['True', 'False']])
    axs[0].set_ylabel('Proportion (%)', labelpad=0, fontdict={'fontsize': LABEL_SIZE})
    axs[0].set_yticks([0, 0.5, 1], [0, 50, 100], fontdict={'fontsize': ANNOT_SIZE})
    axs[0].set_xticks([0, 1], ['Expr. addiction', 'Other genes'], fontdict={'fontsize': ANNOT_SIZE}, rotation=30, ha='right')
    axs[0].set_ylim(0, 1)
    axs[0].spines[['top', 'right']].set_visible(False)
    axs[0].set_title('Oncogenes', fontdict={'fontsize': LABEL_SIZE}, pad=3)
    axs[0].text(0, oncogene_crosstab_pct_df.loc['False', 'True'] / 2, oncogene_crosstab_df.loc['False', 'True'], 
                ha='center', va='center', fontdict={'color': 'white', 'fontsize': ANNOT_SIZE})
    axs[0].text(0, 1 - (oncogene_crosstab_pct_df.loc['True', 'True'] / 2) - LABEL_VERTICAL_ADJUST, oncogene_crosstab_df.loc['True', 'True'], 
                ha='center', va='center', fontdict={'color': 'black', 'fontsize': ANNOT_SIZE})
    axs[0].text(1, oncogene_crosstab_pct_df.loc['False', 'False'] / 2, oncogene_crosstab_df.loc['False', 'False'], 
                ha='center', va='center', fontdict={'color': 'white', 'fontsize': ANNOT_SIZE})
    axs[0].text(1, 1 - (oncogene_crosstab_pct_df.loc['True', 'False'] / 2) - LABEL_VERTICAL_ADJUST, oncogene_crosstab_df.loc['True', 'False'], 
                ha='center', va='center', fontdict={'color': 'black', 'fontsize': ANNOT_SIZE})
    offset_x_labels(fig, axs[0], 25, 0)

    # plot the TFs in the right subplot, with percentages on the axes and numbers in the annotations
    axs[1].bar(['Expr. addiction', 'Other genes'], tf_crosstab_pct_df.loc['False', ['True', 'False']])
    axs[1].bar(['Expr. addiction', 'Other genes'], tf_crosstab_pct_df.loc['True', ['True', 'False']], 
               bottom=tf_crosstab_pct_df.loc['False', ['True', 'False']])
    axs[1].set_xticks([0, 1], ['Expr. addiction', 'Other genes'], fontdict={'fontsize': ANNOT_SIZE}, rotation=30, ha='right')
    axs[1].spines[['top', 'right']].set_visible(False)
    axs[1].set_title('TFs', fontdict={'fontsize': LABEL_SIZE}, pad=3)
    axs[1].text(0, tf_crosstab_pct_df.loc['False', 'True'] / 2, tf_crosstab_df.loc['False', 'True'], 
                ha='center', va='center', fontdict={'color': 'white', 'fontsize': ANNOT_SIZE})
    axs[1].text(0, 1 - (tf_crosstab_pct_df.loc['True', 'True'] / 2) - LABEL_VERTICAL_ADJUST, tf_crosstab_df.loc['True', 'True'], 
                ha='center', va='center', fontdict={'color': 'black', 'fontsize': ANNOT_SIZE})
    axs[1].text(1, tf_crosstab_pct_df.loc['False', 'False'] / 2, tf_crosstab_df.loc['False', 'False'], 
                ha='center', va='center', fontdict={'color': 'white', 'fontsize': ANNOT_SIZE})
    axs[1].text(1, 1 - (tf_crosstab_pct_df.loc['True', 'False'] / 2) - LABEL_VERTICAL_ADJUST, tf_crosstab_df.loc['True', 'False'], 
                ha='center', va='center', fontdict={'color': 'black', 'fontsize': ANNOT_SIZE})
    offset_x_labels(fig, axs[1], 25, 0)
    plt.savefig(os.path.join(FIGURE_DIR, 'Fig_2b_next_gen_expr_addiction_enrichment_barplot.pdf'))

# figure 2c

def plot_expression_addiction_example(dep, x_range=None, y_range=None):
    """
    Plot scatterplot of dependency vs. expression for example expression addictions

    Args:
        dep (str): The name of the gene to plot
        x_range (list[float], optional): x-axis tick values to label. Defaults to None, indicating automatic labels
        y_range (list[float], optional): y-axis tick values to label. Defaults to None, indicating automatic labels
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
    
    # map the gene names to identifiers
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
    next_gen_example_df = construct_example_df(dep_id, dep_id, next_gen_screen_metadata)
    adherent_example_df = construct_example_df(dep_id, dep_id, adherent_lin_matched_screen_metadata)
    
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

    # label axes
    axs[0].set_ylabel('')
    axs[1].set_ylabel('')
    fig.supylabel(f'          {dep} Dependency', fontsize=LABEL_SIZE, x=0.03)
    if y_range is None:
        y_range = [-2.5, -2, -1.5, -1, -0.5, 0, 0.5]
    axs[0].set_yticks(y_range, [], fontdict={'size': TICK_SIZE})
    axs[1].set_yticks(y_range, y_range, fontdict={'size': TICK_SIZE})

    axs[0].set_xlabel('')
    axs[1].set_xlabel(f'{dep} Expression', fontdict={'size': LABEL_SIZE})
    if x_range is None:
        x_range = [0, 2, 4, 6, 8, 10]
    axs[0].set_xticks(x_range, [], fontdict={'size': TICK_SIZE})
    axs[1].set_xticks(x_range, x_range, fontdict={'size': TICK_SIZE})

    # add axes
    if min(x_range) <= 0:
        axs[0].axvline(0, color='black', linestyle='dotted', linewidth=0.5)
        axs[1].axvline(0, color='black', linestyle='dotted', linewidth=0.5)
    axs[0].axhline(0, color='black', linestyle='dotted', linewidth=0.5)
    axs[1].axhline(0, color='black', linestyle='dotted', linewidth=0.5)

    # set common axis limits
    minx = np.min([axs[0].get_xlim()[0], axs[1].get_xlim()[0]])
    maxx = np.max([axs[0].get_xlim()[1], axs[1].get_xlim()[1]])
    miny = np.min([axs[0].get_ylim()[0], axs[1].get_ylim()[0]])
    maxy = np.max([axs[0].get_ylim()[1], axs[1].get_ylim()[1]])
    axs[0].set_xlim(minx, maxx)
    axs[1].set_xlim(minx, maxx)
    axs[0].set_ylim(miny, maxy)
    axs[1].set_ylim(miny, maxy)

    # label each figure with the dataset used
    next_gen_r = next_gen_example_df[['DependencyGeneEffect', 'BiomarkerExpression']].dropna().corr().iloc[0, 1]
    adherent_r = adherent_example_df[['DependencyGeneEffect', 'BiomarkerExpression']].dropna().corr().iloc[0, 1]
    if 'HNRNPH1' in dep:
        axs[0].text(0.05, 0.16, 'NextGen', fontdict={'fontsize': ANNOT_SIZE, 'color': 'black'}, horizontalalignment='left', transform=axs[0].transAxes)
        axs[0].text(0.05, 0.05, 'r = {:.2f}'.format(next_gen_r), 
                    fontdict={'fontsize': ANNOT_SIZE}, ha='left', transform=axs[0].transAxes)
        axs[1].text(0.05, 0.16, 'Traditional', fontdict={'fontsize': ANNOT_SIZE, 'color': 'black'}, horizontalalignment='left', transform=axs[1].transAxes)
        axs[1].text(0.05, 0.05, 'r = {:.2f}'.format(adherent_r), 
                    fontdict={'fontsize': ANNOT_SIZE}, ha='left', transform=axs[1].transAxes)
    elif 'ESR1' in dep:
        axs[0].text(0.95, 0.88, 'NextGen', fontdict={'fontsize': ANNOT_SIZE, 'color': 'black'}, horizontalalignment='right', transform=axs[0].transAxes)
        axs[0].text(0.1, 0.05, 'r = {:.2f}'.format(next_gen_r), 
                    fontdict={'fontsize': ANNOT_SIZE}, ha='left', transform=axs[0].transAxes)
        axs[1].text(0.95, 0.88, 'Traditional', fontdict={'fontsize': ANNOT_SIZE, 'color': 'black'}, horizontalalignment='right', transform=axs[1].transAxes)
        axs[1].text(0.1, 0.05, 'r = {:.2f}'.format(adherent_r), 
                    fontdict={'fontsize': ANNOT_SIZE}, ha='left', transform=axs[1].transAxes)

    plt.savefig(os.path.join(FIGURE_DIR, f'Fig_2c_{dep}_expr_addiction_example.pdf'))
    
    example_summary_df = pd.concat([
        next_gen_example_df.assign(ScreenSet='NextGen'),
        adherent_example_df.assign(ScreenSet='Traditional2D')
    ], axis=0).rename_axis('ScreenID')
    example_summary_df.to_csv(os.path.join(PROCESSED_DIR, f'{dep}_expr_addiction_example_table.csv'))
    
def main():
    """
    Generate all plots related to expressino addictions
    """
    # figure 2a
    plot_expression_addiction_volcano(fdr_threshold=0.005, correlation_threshold=0.3, manual_annotate=MANUAL_ANNOTATE)
    
    # figure 2b
    plot_expression_addiction_enrichment_barplots()

    # figure 2c
    plot_expression_addiction_example('HNRNPH1', x_range=[6, 8, 10], y_range=[-4, -3, -2, -1, 0])
    plot_expression_addiction_example('ESR1', x_range=[0, 2, 4, 6], y_range=[-2, -1, 0, 1])
    

if __name__ == "__main__":
    main()