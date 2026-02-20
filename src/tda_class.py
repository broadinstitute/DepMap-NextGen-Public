import os
import pandas as pd
import numpy as np
from figure_utils import *
from constants import *
from data_utils import *

import seaborn as sns
import matplotlib.pyplot as plt

mpl.style.use(os.path.join(ASSETS_DIR, 'stylesheet.mplstyle'))

def compute_odds_ratio(binary_df, column1, column2, ingroups=[True, True]):
    contingency_table = binary_df.value_counts([column1, column2]).rename('n').reset_index().pivot(index=column1, columns=column2, values='n')
    a = contingency_table.loc[ingroups[0], ingroups[1]]
    b = contingency_table.loc[ingroups[0], not ingroups[1]]
    c = contingency_table.loc[not ingroups[0], ingroups[1]]
    d = contingency_table.loc[not ingroups[0], not ingroups[1]]
    return (a * d) / (b * c)

def run_fishers_exact_from_binary(binary_df, column1, column2, ingroups=[True, True], alternative='two-sided'):    
    contingency_table = binary_df.value_counts([column1, column2]).rename('n').reset_index().pivot(index=column1, columns=column2, values='n')
    a = contingency_table.loc[ingroups[0], ingroups[1]]
    b = contingency_table.loc[ingroups[0], not ingroups[1]]
    c = contingency_table.loc[not ingroups[0], ingroups[1]]
    d = contingency_table.loc[not ingroups[0], not ingroups[1]]
    
    pval = scipy.stats.fisher_exact(table=[[a, b], [c, d]], alternative=alternative).pvalue
    
    return pval

def load_data():
    tda_summary = load_file('tda_dependency_classification.csv', local_dir=PROCESSED_DIR, index_col=0)
    screen_gene_effect = load_file('screen_gene_effect.csv', index_col=0)
    screen_metadata = load_file('screen_metadata.csv', index_col=0)
    geneset_table = load_file('geneset_table_updated.csv')

    gene_functional_class_df = geneset_table[
        geneset_table['Geneset'].isin(['Oncogenes', 'Tumor suppressors', 'Oncology drug targets'])
    ].assign(indicator=True).pivot(index='cds_gene_id', columns='Geneset', values='indicator')

    tda_class_enrichment_df = tda_summary[['NextGenClass']].assign(indicator=True).pivot(
        columns='NextGenClass', values='indicator').fillna(False).join(
        gene_functional_class_df
    ).fillna(False)
    tda_class_enrichment_df = tda_class_enrichment_df[~tda_class_enrichment_df['Unscreened']].drop('Unscreened', axis=1)
    tda_class_enrichment_df['SelOrHighVar'] = tda_class_enrichment_df['Strongly Selective'] | tda_class_enrichment_df['High Variance']
    
    return tda_summary, screen_gene_effect, screen_metadata, gene_functional_class_df, tda_class_enrichment_df

def compute_class_enrichments(tda_class_enrichment_df):
    tda_classes = [
        # 'Growth Suppressor',
        'Non-dependency',
        'Weakly Selective',
        'Strongly Selective',
        'High Variance',
        'Pan-dependency'
    ]
    functional_classes = [
        'Oncogenes',
        'Tumor suppressors',
        'Oncology drug targets'
    ]
    
    tda_enrichment_matrix = np.zeros((len(tda_classes), len(functional_classes)))
    tda_class_fisher_p_matrix = np.ones_like(tda_enrichment_matrix)
    for i, c1 in enumerate(tda_classes):
        for j, c2 in enumerate(functional_classes):
            tda_enrichment_matrix[i, j] = compute_odds_ratio(tda_class_enrichment_df, c1, c2, ingroups=[True, True])
            tda_class_fisher_p_matrix[i, j] = run_fishers_exact_from_binary(tda_class_enrichment_df, c1, c2, ingroups=[True, True])

    tda_enrichment_matrix = pd.DataFrame(np.log2(tda_enrichment_matrix), index=tda_classes, columns=functional_classes)
    tda_class_fisher_p_matrix = pd.DataFrame(tda_class_fisher_p_matrix, index=tda_classes, columns=functional_classes)
    tda_class_fisher_fdr_matrix = tda_class_fisher_p_matrix.melt(ignore_index=False).reset_index().rename({'index': 'tda_class', 'variable': 'func_class', 'value': 'p'}, axis=1).assign(
        FDR=lambda x: scipy.stats.false_discovery_control(x['p'])
    ).pivot(index='tda_class', columns='func_class', values='FDR').reindex_like(tda_class_fisher_p_matrix)
    tda_class_significance_strings = tda_class_fisher_fdr_matrix.applymap(format_significance).reindex_like(tda_class_fisher_p_matrix).replace({'n.s.': ''})
    
    formatted_enrichment_df = tda_enrichment_matrix.melt(ignore_index=False).reset_index().rename({'index': 'NextGenClass', 'variable': 'FunctionalClass', 'value': 'OddsRatio'}, axis=1).merge(
        tda_class_fisher_p_matrix.melt(ignore_index=False).reset_index().rename({'index': 'NextGenClass', 'variable': 'FunctionalClass', 'value': 'FisherP'}, axis=1).assign(
            FisherFDR=lambda x: scipy.stats.false_discovery_control(x['FisherP'])
        )
    )
    
    return tda_enrichment_matrix, tda_class_fisher_p_matrix, tda_class_fisher_fdr_matrix, tda_class_significance_strings, formatted_enrichment_df

def print_selected_enrichments(tda_class_enrichment_df, tda_enrichment_matrix, tda_class_fisher_fdr_matrix):
    print('Oncology drug targets in strongly selective + high variance:')
    print('OR:', compute_odds_ratio(tda_class_enrichment_df, 'SelOrHighVar', 'Oncology drug targets', ingroups=[True, True]))
    print("Fisher's p:", run_fishers_exact_from_binary(tda_class_enrichment_df, 'SelOrHighVar', 'Oncology drug targets', ingroups=[True, True]))
    
    print('Oncogenes in strongly selective + high variance:')
    print('OR:', compute_odds_ratio(tda_class_enrichment_df, 'SelOrHighVar', 'Oncogenes', ingroups=[True, True]))
    print("Fisher's p:", run_fishers_exact_from_binary(tda_class_enrichment_df, 'SelOrHighVar', 'Oncogenes', ingroups=[True, True]))

    print('Tumor suppressor genes in strongly selective + high variance:')
    print('OR:', compute_odds_ratio(tda_class_enrichment_df, 'SelOrHighVar', 'Tumor suppressors', ingroups=[True, True]))
    print("Fisher's p:", run_fishers_exact_from_binary(tda_class_enrichment_df, 'SelOrHighVar', 'Tumor suppressors', ingroups=[True, True]))

# supplementary

def plot_tda_class_distributions_with_enrichment(screen_metadata, screen_gene_effect, tda_summary, tda_enrichment_matrix, tda_class_significance_strings):
    
    next_gen_screens_list = screen_metadata[
        screen_metadata['ScreenType'].isin(['3DO', '2DN', '3DN'])
    ].index.tolist()
    next_gen_screens_list = list(set(next_gen_screens_list) & set(screen_gene_effect.index))
    
    ##### NOTE: we may need to manually select representative genes from each class which have visually distinct dependency profiles #####
    representative_genes = ['RHCE (6006)', 'SLC35C1 (55343)', 'ATP6V1G2 (534)', 'ITGB1 (3688)', 'CUL1 (8454)']
    
    dependency_class_summary_df = screen_gene_effect.loc[next_gen_screens_list, representative_genes].melt(
        var_name='cds_gene_id', value_name='gene_effect', ignore_index=False
    ).reset_index().merge(tda_summary['NextGenClass'], left_on='cds_gene_id', right_index=True)
    dependency_class_size = tda_summary['NextGenClass'].value_counts().to_dict()
    
    tda_classes = tda_enrichment_matrix.index.tolist()
    ncat = len(tda_classes)
    
    fig, ax_map = prepare_gridspec(
        ax_map=pd.DataFrame([[True]]),
        figsize=(55 * mm, 45 * mm),
        height_ratios=[1], wspace=0.7,
        inner_gs_dict={(0, 0): {'nrows': ncat, 'ncols':1, 'hspace': 0}}
    )
    plt.subplots_adjust(left=0.35, bottom=0.28, right=0.95, top=0.85)

    # plot the gene effect distributions as densities
    xlims = np.zeros((ncat, 2))
    for i, ms in enumerate(tda_classes):
        # subset = dependency_class_summary_df[dependency_class_summary_df['OrganoidsNeurospheresClass'] == ms]
        subset = dependency_class_summary_df[dependency_class_summary_df['NextGenClass'] == ms]
        sns.kdeplot(subset, x='gene_effect', color=dependency_class_palette[ms], fill=True, linewidth=0.5, ax=ax_map[0][0][(i, 0)])
        ax_map[0][0][(i, 0)].axhline(0, c='black')
        xlims[i, :] = ax_map[0][0][(i, 0)].get_xlim()

    minx = -3.25 #np.min(xlims) 
    maxx = 1.25 #np.max(xlims)

    # reformat all plots and remove extra legends
    for i, ms in enumerate(tda_classes):
        ax_map[0][0][(i, 0)].set_xlim(minx, maxx)
        ax_map[0][0][(i, 0)].set_ylim(0, 1.1 * ax_map[0][0][(i, 0)].get_ylim()[1])

        ax_map[0][0][(i, 0)].set(yticks=[])
        ax_map[0][0][(i, 0)].set_ylabel(ms + f'\nn = {dependency_class_size[ms]:,}', rotation=0, color='black', labelpad=4, 
                                horizontalalignment='right', verticalalignment='center', size=TICK_SIZE)
        ax_map[0][0][(i, 0)].set_xlabel('')

        ax_map[0][0][(i, 0)].axvline(0, color='black', linestyle='dotted', alpha=0.5, zorder=-1, linewidth=0.5)

        if i == ncat - 1:
            ax_map[0][0][(i, 0)].set_xlabel('Mean CRISPR dependency\n(Chronos)', size=LABEL_SIZE)
            ax_map[0][0][(i, 0)].set_xticks([-3, -2, -1, 0, 1])
            ax_map[0][0][(i, 0)].set_xticklabels(['-3', '-2', '-1', '0', '1'], rotation=0, color='black')

        ax_map[0][0][(i, 0)].spines[['bottom', 'top']].set_linewidth(0)
        ax_map[0][0][(i, 0)].spines[['left', 'right']].set_linewidth(0.5)
    ax_map[0][0][(0, 0)].spines[['top']].set_linewidth(0.5)
    ax_map[0][0][(len(tda_classes)-1, 0)].spines[['bottom']].set_linewidth(0)

    plt.savefig(os.path.join(FIGURE_DIR, 'tda_dependency_class_breakdown.pdf'))
    return dependency_class_summary_df.rename({'index': 'ScreenID', 'cds_gene_id': 'Gene', 'gene_effect': 'GeneEffect'}, axis=1)
    
        
def main():
    tda_summary, screen_gene_effect, screen_metadata, gene_functional_class_df, tda_class_enrichment_df = load_data()
    
    tda_enrichment_matrix, tda_class_fisher_p_matrix, tda_class_fisher_fdr_matrix, tda_class_significance_strings, formatted_enrichment_df = compute_class_enrichments(tda_class_enrichment_df)
    
    print_selected_enrichments(tda_class_enrichment_df, tda_enrichment_matrix, tda_class_fisher_fdr_matrix)

    # supplementary
    representative_gene_effects_df = plot_tda_class_distributions_with_enrichment(screen_metadata, screen_gene_effect, tda_summary, tda_enrichment_matrix, tda_class_significance_strings)
    
    formatted_enrichment_df.to_csv(os.path.join(PROCESSED_DIR, 'dependency_class_enrichment.csv'), index=False)
    representative_gene_effects_df.to_csv(os.path.join(PROCESSED_DIR, 'tda_class_representative_gene_effects.csv'), index=False)
    

if __name__ == "__main__":
    main()