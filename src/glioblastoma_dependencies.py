import os
import pandas as pd
import numpy as np
from figure_utils import *
from constants import *
from data_utils import *

import seaborn as sns
import matplotlib.pyplot as plt

mpl.style.use(os.path.join(ASSETS_DIR, 'stylesheet.mplstyle'))

MANUAL_ANNOTATE = False

# figure ED 6f

def plot_cdk6_dependency_expression_volcano(fdr_threshold=0.005, correlation_threshold=0.3, highlight_color='#C54B8C', manual_annotate=MANUAL_ANNOTATE):
    cdk6_expr_volcano_df = load_file('gbm_cdk6_dependency_expression_tests.csv', local_dir=PROCESSED_DIR, index_col=0)
    cdk6_expr_volcano_df['p_transform'] = -np.log10(cdk6_expr_volcano_df['p-value'])
    cdk6_expr_volcano_df['label'] = cdk6_expr_volcano_df.index.str.split(' ').str[0]
    
    xlabel='Pearson Correlation'
    ylabel='p_transform'
    
    # make the plot
    plt.figure(figsize=(45 * mm, 45 * mm))
    plt.subplots_adjust(left=0.22, bottom=0.22)
    
    sns.scatterplot(
        cdk6_expr_volcano_df.sort_values('FDR', ascending=False),
        x=xlabel, y=ylabel,
        hue='Significant', palette={True: highlight_color, False: 'tab:gray'},
        s=4,
        alpha=0.6, edgecolor='black', linewidth=0.25, legend=False
    )
    plt.xlabel('Pearson correlation ($r$)', fontdict={'size': LABEL_SIZE})
    plt.ylabel('$-\log_{10}(P)$', fontdict={'size': LABEL_SIZE})
    # plt.title('Expression Addictions', fontdict={'size': TITLE_SIZE})
    
    # add a significance line
    p_fdr_approx = approximate_fdr_in_p_units(cdk6_expr_volcano_df, fdr_threshold=fdr_threshold, pcolumn='p-value', fdr_column='FDR')
    plt.axhline(-np.log10(p_fdr_approx), linestyle='dashed', color='tab:red', alpha=0.5, linewidth=0.5)
    plt.axvline(-correlation_threshold, linestyle='dashed', color='tab:red', alpha=0.5, linewidth=0.5)
    plt.axvline(correlation_threshold, linestyle='dashed', color='tab:red', alpha=0.5, linewidth=0.5)
    
    
    # annotate the top genes
    annots = cdk6_expr_volcano_df.set_index('label')
    
    if manual_annotate:
        # top genes
        manually_annotate(plt.gca(), 'NOS2', annots, 
                          xlabel, ylabel, (0.05, 0), ha='left', va='center', color=highlight_color,
                          arrowprops=dict(arrowstyle='-', color='black', alpha=1, linewidth=0.25, shrinkA=0, shrinkB=1),
                          bbox_kws={'facecolor': 'white', 'edgecolor': 'white', 'boxstyle':'round,pad=0.05'})
        manually_annotate(plt.gca(), 'SLC1A5', annots, 
                          xlabel, ylabel, (-0.05, 0.5), ha='right', va='center', color=highlight_color,
                          arrowprops=dict(arrowstyle='-', color='black', alpha=1, linewidth=0.25, shrinkA=0, shrinkB=1),
                          bbox_kws={'facecolor': 'white', 'edgecolor': 'white', 'boxstyle':'round,pad=0.05'})
        # selected genes
        manually_annotate(plt.gca(), 'CDKN2A', annots, 
                          xlabel, ylabel, (-0.075, -0.21), ha='right', va='center', color=highlight_color,
                          arrowprops=dict(arrowstyle='-', color='black', alpha=1, linewidth=0.25, shrinkA=0, shrinkB=1),
                          bbox_kws={'facecolor': 'white', 'edgecolor': 'white', 'boxstyle':'round,pad=0.05'})
        manually_annotate(plt.gca(), 'CDK6', annots, 
                          xlabel, ylabel, (0.06, 0.7), ha='left', va='center', color=highlight_color,
                          arrowprops=dict(arrowstyle='-', color='black', alpha=1, linewidth=0.25, shrinkA=0, shrinkB=1),
                          bbox_kws={'facecolor': 'white', 'edgecolor': 'white', 'boxstyle':'round,pad=0.05'})
        manually_annotate(plt.gca(), 'MET', annots, 
                          xlabel, ylabel, (-0.07, 0.8), ha='right', va='center', color=highlight_color,
                          arrowprops=dict(arrowstyle='-', color='black', alpha=1, linewidth=0.25, shrinkA=0, shrinkB=1),
                          bbox_kws={'facecolor': 'white', 'edgecolor': 'white', 'boxstyle':'round,pad=0.05'})
        manually_annotate(plt.gca(), 'EGFR', annots, 
                          xlabel, ylabel, (-0.25, 3), ha='left', va='center', color='tab:gray',
                          arrowprops=dict(arrowstyle='-', color='black', alpha=1, linewidth=0.25, shrinkA=0, shrinkB=1),
                          bbox_kws={'facecolor': 'white', 'edgecolor': 'white', 'boxstyle':'round,pad=0.05'})
        manually_annotate(plt.gca(), 'PTEN', annots, 
                          xlabel, ylabel, (-0.05, 1.8), ha='left', va='center', color='tab:gray',
                          arrowprops=dict(arrowstyle='-', color='black', alpha=1, linewidth=0.25, shrinkA=0, shrinkB=1),
                          bbox_kws={'facecolor': 'white', 'edgecolor': 'white', 'boxstyle':'round,pad=0.05'})
    else:
        texts = [plt.gca().text(r[xlabel], r[ylabel], i, ha='center', va='center', fontsize=ANNOT_SIZE) for i,r in annots.head(10).iterrows()]
    
    plt.savefig(os.path.join(FIGURE_DIR, f'gbm_cdk6_dependency_expression_volcano.pdf'))

# figure ED 6g

def plot_cdk6_dependency_mutation_volcano(fdr_threshold=0.05, mean_difference_threshold=0.3, manual_annotate=MANUAL_ANNOTATE):
    cdk6_mutation_volcano_df = load_file('gbm_cdk6_dependency_mutation_tests.csv', local_dir=PROCESSED_DIR)
    
    bm_type_map = {'PTEN': '\ndamaging', 'TP53': '\ndamaging', 'NF1': '\ndamaging', 'PIK3CA': '\nhotspot', 'TERT': 'promoter\nhotspot', 'BRAF': '\nhotspot'}
    
    cdk6_mutation_volcano_df['-log10(MWU p)'] = -np.log10(cdk6_mutation_volcano_df['MWU p'])
    cdk6_mutation_volcano_df['label'] = cdk6_mutation_volcano_df['Biomarker'].str.split(' ').str[0].apply(lambda x: x + ' ' + bm_type_map[x] if x in bm_type_map else x)
    
    xlabel='Mean Difference'
    ylabel='-log10(MWU p)'
    
    plt.figure(figsize=(75 * mm, 50 * mm))
    plt.subplots_adjust(left=0.12, bottom=0.22, right=0.65, top=0.84)

    sns.scatterplot(
        cdk6_mutation_volcano_df,
        x=xlabel, y=ylabel,
        hue='Significant', palette={True: '#C54B8C', False: 'tab:gray'},
        s=8,
        alpha=0.6, edgecolor='black', legend=False
    )
    
    p_fdr_approx = approximate_fdr_in_p_units(cdk6_mutation_volcano_df, fdr_threshold=fdr_threshold, pcolumn='MWU p', fdr_column='FDR')
    # plt.axhline(-np.log10(p_fdr_approx), linestyle='dashed', color='tab:red', alpha=0.5, linewidth=0.5)
    plt.axvline(mean_difference_threshold, linestyle='dashed', color='tab:red', alpha=0.5, linewidth=0.5)
    plt.axvline(-mean_difference_threshold, linestyle='dashed', color='tab:red', alpha=0.5, linewidth=0.5)

    plt.ylabel('$-\log_{10}(P)$', fontdict={'size': LABEL_SIZE})
    plt.xlabel('Mean Dependency Difference\n(Altered - Unaltered)', fontdict={'size': LABEL_SIZE})
    
    annots = cdk6_mutation_volcano_df.set_index('label')
    
    if manual_annotate:
        manually_annotate(plt.gca(), 'PTEN \ndamaging', annots, 
                          xlabel, ylabel, (0, 0.04), ha='center', va='center',
                          arrowprops=dict(arrowstyle='-', color='black', alpha=0, linewidth=0.25, shrinkA=0, shrinkB=1),
                          bbox_kws={'facecolor': 'white', 'edgecolor': 'white', 'boxstyle':'round,pad=0.05'})
        manually_annotate(plt.gca(), 'TP53 \ndamaging', annots, 
                          xlabel, ylabel, (-0.02, -0.03), ha='right', va='center',
                          arrowprops=dict(arrowstyle='-', color='black', alpha=0, linewidth=0.25, shrinkA=0, shrinkB=1),
                          bbox_kws={'facecolor': 'white', 'edgecolor': 'white', 'boxstyle':'round,pad=0.05'})
        manually_annotate(plt.gca(), 'NF1 \ndamaging', annots, 
                          xlabel, ylabel, (-0.02, 0), ha='right', va='center',
                          arrowprops=dict(arrowstyle='-', color='black', alpha=0, linewidth=0.25, shrinkA=0, shrinkB=1),
                          bbox_kws={'facecolor': 'white', 'edgecolor': 'white', 'boxstyle':'round,pad=0.05'})
        manually_annotate(plt.gca(), 'PIK3CA \nhotspot', annots, 
                          xlabel, ylabel, (-0.02, 0), ha='right', va='center',
                          arrowprops=dict(arrowstyle='-', color='black', alpha=0, linewidth=0.25, shrinkA=0, shrinkB=1),
                          bbox_kws={'facecolor': 'white', 'edgecolor': 'white', 'boxstyle':'round,pad=0.05'})
        manually_annotate(plt.gca(), 'TERT promoter\nhotspot', annots, 
                          xlabel, ylabel, (-0.02, 0), ha='right', va='center',
                          arrowprops=dict(arrowstyle='-', color='black', alpha=0, linewidth=0.25, shrinkA=0, shrinkB=1),
                          bbox_kws={'facecolor': 'white', 'edgecolor': 'white', 'boxstyle':'round,pad=0.05'})
        manually_annotate(plt.gca(), 'BRAF \nhotspot', annots, 
                          xlabel, ylabel, (-0.02, 0), ha='right', va='center',
                          arrowprops=dict(arrowstyle='-', color='black', alpha=0, linewidth=0.25, shrinkA=0, shrinkB=1),
                          bbox_kws={'facecolor': 'white', 'edgecolor': 'white', 'boxstyle':'round,pad=0.05'})
    else:
        texts = [plt.gca().text(r[xlabel], r[ylabel], i, ha='center', va='center', fontsize=ANNOT_SIZE) for i,r in annots.head(10).iterrows()]


    plt.savefig(os.path.join(FIGURE_DIR, 'gbm_cdk6_dependency_mutation_volcano.pdf'))
    

def main():
    # figure ED6f
    plot_cdk6_dependency_expression_volcano(fdr_threshold=0.005, correlation_threshold=0.3, manual_annotate=MANUAL_ANNOTATE)
    
    # figure ED6f
    plot_cdk6_dependency_mutation_volcano(fdr_threshold=0.05, mean_difference_threshold=0.3, manual_annotate=MANUAL_ANNOTATE)
    

if __name__ == "__main__":
    main()
