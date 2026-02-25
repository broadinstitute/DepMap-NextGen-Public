import os
import pandas as pd
import numpy as np
from figure_utils import *
from constants import *
from data_utils import *
from gene_utils import *

import seaborn as sns
import matplotlib.pyplot as plt

mpl.style.use(os.path.join(ASSETS_DIR, 'stylesheet.mplstyle'))

MANUAL_ANNOTATE = True

# reported stats

def print_combined_enrichment(top_n_each=200):
    """
    Print enrichment significance for the cluster of adhesion gene sets.

    Args:
        top_n_each (int, optional): Number of genes to consider for overrepresentation analysis. Defaults to 200.
    """
    org_2d_diff_geneset_df = load_file('organoids_vs_2d_geneset_enrichment_table.csv', local_dir=PROCESSED_DIR)
    geneset_clusters = org_2d_diff_geneset_df.set_index('Geneset')['Cluster']
    org_2d_diff_df = load_file('organoids_vs_2d_dependency_volcano_table.csv', local_dir=PROCESSED_DIR, index_col=0)
    geneset_table = load_file('geneset_table_updated.csv')
    combined_geneset_dict = geneset_table[geneset_table['Source'].isin(['Hallmark', 'KEGG'])].groupby('Geneset').apply(lambda x: list(x['GeneSymbol'])).to_dict()

    representative_adhesion_set = 'KEGG_FOCAL_ADHESION'
    
    adhesion_genesets = geneset_clusters.loc[lambda x: x == geneset_clusters.loc[representative_adhesion_set]].index.tolist()
    combined_adhesion_geneset = list(set.union(*[set(combined_geneset_dict[gs]) for gs in adhesion_genesets]))
    gene_universe = org_2d_diff_df['GeneSymbol'].tolist()
    effective_adhesion_geneset = list(set(combined_adhesion_geneset) & set(gene_universe))

    combined_adhesion_test = run_single_hypergeometric(
        org_2d_diff_df.sort_values('MeanDifference')['GeneSymbol'].head(top_n_each).tolist() + org_2d_diff_df.sort_values('MeanDifference')['GeneSymbol'].tail(top_n_each).tolist(), 
        combined_adhesion_geneset, 
        org_2d_diff_df['GeneSymbol'].tolist()
    )
    print('Adhesion cluster combined hypergeometric: ')
    print('OR:', combined_adhesion_test.odds_ratio)
    print('p-value:', combined_adhesion_test.pval)

# figure ED10e

def plot_geneset_similarity_matrix():
    """
    Plot the similarity between genesets as a heatmap
    """
    sig_geneset_similarity = load_file('organoids_vs_2d_top_geneset_similarity_matrix.csv', local_dir=PROCESSED_DIR, index_col=0)
    
    plt.figure(figsize=(90 * mm, 40 * mm))
    plt.subplots_adjust(left=0.5, bottom=0.1, top=0.9, right=0.92)
    sns.heatmap(
        sig_geneset_similarity, cmap='Reds', vmin=0, vmax=1,
        xticklabels=False,
        cbar_kws={'label': 'Overlap coefficient'}
    )
    plt.savefig(os.path.join(FIGURE_DIR, 'Fig_ED12e_organoids_vs_2d_geneset_similarity.pdf'))
    
# figure 5a

def plot_geneset_stemplot():
    """
    Plot the enrichment of genesets within the tails of the distribution of dependency differences between
    Organoid and Traditional 2D models, colored by geneset cluster
    """
    org_2d_diff_geneset_df = load_file('organoids_vs_2d_geneset_enrichment_table.csv', local_dir=PROCESSED_DIR)
    geneset_clusters = org_2d_diff_geneset_df.set_index('Geneset')['Cluster']

    org_2d_diff_geneset_df = org_2d_diff_geneset_df.sort_values('FDR')
    org_2d_diff_geneset_df['clean_term'] = org_2d_diff_geneset_df['Geneset'].apply(
        lambda x: clean_geneset_name(x, nth=10, remove_prefix=False)
    )
    org_2d_diff_geneset_df['-log10(FDR)'] = -np.log10(org_2d_diff_geneset_df['FDR'])
    org_2d_diff_geneset_df['Cluster'] = org_2d_diff_geneset_df['Cluster'].astype(str)
    
    
    plt.figure(figsize=(90 * mm, 45 * mm))
    plt.subplots_adjust(left=0.5, bottom=0.15, top=0.85, right=0.95)

    # plot the stems and the heads separately
    sns.barplot(
        org_2d_diff_geneset_df.sort_values('FDR'),
        x='OddsRatio',
        y='clean_term',
        hue='Cluster',
        palette=cluster_cmap,
        dodge=False,
        legend=False,
        width=0.1,
        alpha=0.7
    )
    sns.scatterplot(
        org_2d_diff_geneset_df.sort_values('FDR'),
        x='OddsRatio',
        y='clean_term',
        hue='Cluster',
        palette=cluster_cmap,
        size='-log10(FDR)',
        sizes=(8, 24),
        edgecolor='black'
    )

    # trim the legend to use 3 sizes instead of 5
    handles = plt.gca().get_legend_handles_labels()[0][len(geneset_clusters.unique())+1:]
    smallest_size, smallest_label = handles[1].get_markersize(), np.round(float(handles[1].get_label()), 2)
    largest_size, largest_label = handles[-1].get_markersize(), np.round(float(handles[-1].get_label()), 2)
    median_size, median_label = (smallest_size + largest_size) / 2, (smallest_label + largest_label) / 2

    new_handles = [
        Patch(facecolor='white', edgecolor='white', label='-log10(FDR)'),
        Line2D([], [], color='black', marker = 'o', linestyle='None', markeredgewidth=0.5,
                          markersize=smallest_size, label = str(np.round(smallest_label, 2))),
        Line2D([], [], color='black', marker = 'o', linestyle='None', markeredgewidth=0.5,
                          markersize=median_size, label = str(np.round(median_label, 2))),
        Line2D([], [], color='black', marker = 'o', linestyle='None', markeredgewidth=0.5,
                          markersize=largest_size, label = str(np.round(largest_label, 2)))
    ]
    plt.legend(handles=new_handles, prop={'size': ANNOT_SIZE})
    plt.title('Highly Differing Genesets', fontdict={'size': TITLE_SIZE})
    plt.ylabel('')
    plt.yticks(range(len(org_2d_diff_geneset_df)), org_2d_diff_geneset_df.sort_values('FDR')['clean_term'].tolist(), fontdict={'size': TICK_SIZE})
    plt.xlabel('Odds Ratio', fontdict={'size': LABEL_SIZE})
    plt.xticks(plt.gca().get_xticks()[:-1], labels=plt.gca().get_xticklabels()[:-1], fontdict={'size': TICK_SIZE})
    plt.savefig(os.path.join(FIGURE_DIR, 'Fig_5a_organoid_vs_2d_geneset_enrichment_stemplot.pdf'))
    
# figure 5b
    
def plot_global_volcano_plot(manual_annotate=MANUAL_ANNOTATE):
    """
    Plot the mean difference in dependency between Organoid and Traditional adherent models, colored by geneset cluster

    Args:
        manual_annotate (boolean, optional): Flag for whether to annotate genes manually. Defaults to MANUAL_ANNOTATE.
    """
    org_2d_volcano_df = load_file('organoids_vs_2d_dependency_volcano_table.csv', local_dir=PROCESSED_DIR, index_col=0)
    org_2d_volcano_df['-log10(MWU p)'] = -np.log10(org_2d_volcano_df['MannWhitneyPValue'])
    
    # make the volcano plot
    plt.figure(figsize=(90 * mm, 45 * mm))
    plt.subplots_adjust(bottom=0.22)
    
    xlabel='MeanDifference'
    ylabel='-log10(MWU p)'

    sns.scatterplot(
        org_2d_volcano_df.sort_values(['Cluster', 'MannWhitneyFDR'], ascending=[False, True]),
        x=xlabel,
        y=ylabel,
        hue='Cluster',
        palette=cluster_cmap,
        s=4,
        alpha=0.7,
        legend=False
    )

    # annotate some genes
    annots = org_2d_volcano_df[
        org_2d_volcano_df['GeneSymbol'].isin(['MSMO1', 'SQLE', 'ITGB1', 'ITGAV', 'CDK4', 'GPX4', 'TFRC', 'NSDHL', 'IDI1', 'HSD17B7'])
    ].set_index('GeneSymbol')
    if manual_annotate:
        manually_annotate(plt.gca(), 'MSMO1', annots, 
                          xlabel, ylabel, (0.05, 0), ha='left', va='center',
                          arrowprops=dict(arrowstyle='-', color='black', alpha=1, linewidth=0.25, shrinkA=0, shrinkB=1),
                          color=cluster_cmap[annots.loc['MSMO1', 'Cluster']], weight='bold')
        manually_annotate(plt.gca(), 'SQLE', annots, 
                          xlabel, ylabel, (0.05, 0), ha='left', va='center',
                          arrowprops=dict(arrowstyle='-', color='black', alpha=1, linewidth=0.25, shrinkA=0, shrinkB=1),
                          color=cluster_cmap[annots.loc['SQLE', 'Cluster']], weight='bold')
        manually_annotate(plt.gca(), 'ITGB1', annots, 
                          xlabel, ylabel, (-0.2, 0.1), ha='right', va='center',
                          arrowprops=dict(arrowstyle='-', color='black', alpha=1, linewidth=0.25, shrinkA=0, shrinkB=1),
                          color=cluster_cmap[annots.loc['ITGB1', 'Cluster']], weight='bold')
        manually_annotate(plt.gca(), 'ITGAV', annots, 
                          xlabel, ylabel, (-0.08, 0), ha='left', va='center',
                          arrowprops=dict(arrowstyle='-', color='black', alpha=1, linewidth=0.25, shrinkA=0, shrinkB=1),
                          color=cluster_cmap[annots.loc['ITGAV', 'Cluster']], weight='bold')
        manually_annotate(plt.gca(), 'CDK4', annots, 
                          xlabel, ylabel, (0.07, 3), ha='left', va='center',
                          arrowprops=dict(arrowstyle='-', color='black', alpha=1, linewidth=0.25, shrinkA=0, shrinkB=1),
                          color=cluster_cmap[annots.loc['CDK4', 'Cluster']], weight='normal')
        manually_annotate(plt.gca(), 'GPX4', annots, 
                          xlabel, ylabel, (0.05, 0), ha='left', va='center',
                          arrowprops=dict(arrowstyle='-', color='black', alpha=1, linewidth=0.25, shrinkA=0, shrinkB=1),
                          color=cluster_cmap[annots.loc['GPX4', 'Cluster']], weight='normal')
        manually_annotate(plt.gca(), 'TFRC', annots, 
                          xlabel, ylabel, (0.05, 0), ha='left', va='center',
                          arrowprops=dict(arrowstyle='-', color='black', alpha=1, linewidth=0.25, shrinkA=0, shrinkB=1),
                          color=cluster_cmap[annots.loc['TFRC', 'Cluster']], weight='normal')
        manually_annotate(plt.gca(), 'NSDHL', annots, 
                          xlabel, ylabel, (0.05, 0), ha='left', va='center',
                          arrowprops=dict(arrowstyle='-', color='black', alpha=1, linewidth=0.25, shrinkA=0, shrinkB=1),
                          color=cluster_cmap[annots.loc['NSDHL', 'Cluster']], weight='normal')
        manually_annotate(plt.gca(), 'HSD17B7', annots, 
                          xlabel, ylabel, (0.05, -3), ha='left', va='center',
                          arrowprops=dict(arrowstyle='-', color='black', alpha=1, linewidth=0.25, shrinkA=0, shrinkB=1),
                          color=cluster_cmap[annots.loc['HSD17B7', 'Cluster']], weight='normal')
    else:
        texts = [plt.gca().text(r[xlabel], r[ylabel], i, ha='center', va='center', fontsize=ANNOT_SIZE) for i,r in annots.iterrows()]
    
    plt.title('Global Dependency Differences', fontdict={'size': TITLE_SIZE})
    plt.xlabel('Mean Difference\n(Organoid - 2D)', fontdict={'size': TITLE_SIZE})
    plt.ylabel('-log10(Mann-Whitney U p-value)', fontdict={'size': LABEL_SIZE})
    plt.xticks(plt.gca().get_xticks()[1:-1], labels=plt.gca().get_xticklabels()[1:-1], fontdict={'size': TICK_SIZE})
    plt.yticks(plt.gca().get_yticks()[1:-1], labels=plt.gca().get_yticklabels()[1:-1], fontdict={'size': TICK_SIZE})

    handles = []
    handles.append(Patch(facecolor='white', edgecolor='white', label='Geneset cluster'))
    for gs in ['Lipid Metabolism', 'mTORC1', 'Adhesion/Cytoskeleton', 'OxPhos', 'Cell Cycle', 'SLE', 'UV']:
            handle = Line2D([], [], color=cluster_cmap[gs], marker = 'o', linestyle='None', markersize=4, markeredgewidth=0, label = gs)
            handles.append(handle)    
    plt.legend(handles = handles, ncol=2, handletextpad=0, columnspacing=0.5, loc='upper right', prop={'size': ANNOT_SIZE})

    plt.savefig(os.path.join(FIGURE_DIR, 'Fig_5b_organoid_vs_2d_dependency_volcano.pdf'))
    

def main():
    """
    Generate all plots related to global dependency differences between organoid and adherent models
    """
    print_combined_enrichment()
    
    # figure ED12e
    plot_geneset_similarity_matrix()
    
    # figure 5a
    plot_geneset_stemplot()
    
    # figure 5b
    plot_global_volcano_plot(manual_annotate=MANUAL_ANNOTATE)

    
if __name__ == "__main__":
    main()   
    