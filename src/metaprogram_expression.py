import os
import pandas as pd
import numpy as np
from figure_utils import *
from constants import *
from data_utils import *

import seaborn as sns
import matplotlib.pyplot as plt

mpl.style.use(os.path.join(ASSETS_DIR, 'stylesheet.mplstyle'))

# figure 3a

def plot_metaprogram_bubble():
    """
    Plot the average metaprogram expression per lineage as a bubble (dot) plot.
    """
    # load data and annotations
    mmp_comparison_tests = load_file('next_gen_vs_2d_metaprogram_expression_tests.csv', local_dir=PROCESSED_DIR).replace(lineage_replacement)
    geneset_table = load_file('geneset_table_updated.csv')
    mmp_comparison_tests['Metaprogram'] = mmp_comparison_tests['Metaprogram'].str.strip()
    mmps = geneset_table[(geneset_table['Source'] == 'Gavish')]['Geneset'].str.strip().unique().tolist()
    bubble_lineage_order = ['CNS/Brain', 'All Organoid Lineages', 'Esophagus/Stomach', 'Pancreas', 'Colorectal', 'Breast', 'Ovary/Fallopian Tube', 'Prostate']
    
    # plot dot plot
    axs, cbar = bubble_plot(
        mmp_comparison_tests.assign(fdr_transform=lambda x: -np.log10(x['MannWhitneyFDR'])).rename({'fdr_transform': '-log10(FDR)'}, axis=1),
        x_label='Metaprogram',
        y_label='OncotreeLineage',
        xorder=mmps, 
        yorder=bubble_lineage_order[::-1],
        hue_label='MeanDifference', hue_min=-3, hue_center=0, hue_max=3,
        size_label='-log10(FDR)', size_min=0, size_max=5, sizes=(2, 24),
        linewidth=0.5, edgecolor='black',
        figsize=(180 * mm, 45 * mm), legend_loc=(1, -1.3),
        width_ratios=[59, 1], gridspec_kws={'wspace': 0.15},
    )

    plt.subplots_adjust(left=0.15, bottom=0.5, right=0.9)

    axs[0].set_xlabel('')
    axs[0].set_ylabel('')

    cbar.ax.yaxis.set_label_position('left')
    cbar.ax.set_ylabel('Mean Difference\n(NextGen - Traditional)', rotation=90, labelpad=5, fontsize=LABEL_SIZE)
    cbar.ax.set_yticks([-2, 0, 2], labels=[-2, 0, 2])
    cbar.outline.set_linewidth(0)

    axs[0].set_yticks([i for i in range(len(bubble_lineage_order))], labels=bubble_lineage_order[::-1], fontdict={'size': TICK_SIZE})
    axs[0].set_xticks([i for i in range(len(mmps))], labels=mmps, fontdict={'size': TICK_SIZE})
    axs[0].set_title('NextGen vs. Traditional Metaprogram Expression', fontdict={'size': TITLE_SIZE})

    h, l = axs[0].get_legend_handles_labels()
    axs[0].get_legend().remove()
    axs[0].legend(handles = h[6:7] + h[7::2], loc='lower left', bbox_to_anchor=(1, -0.85), prop={'size': ANNOT_SIZE})

    plt.savefig(os.path.join(FIGURE_DIR, 'Fig_3a_next_gen_vs_2d_metaprogram_expr_diff_bubble.pdf'))
    
# figure 4a

def plot_metaprogram_volcano(mean_difference_threshold=0.5, p_threshold=1e-15):
    """
    Plot the difference in metaprogram expression means between NextGen and Traditional adherent models.

    Args:
        mean_difference_threshold (float, optional): Minimum difference in mean expression score. Defaults to 0.5.
        p_threshold (float, optional): Minimum significance value by Mann-Whitney U / Wilcoxon rank-sum test. Defaults to 1e-15.
    """
    # load data
    mmp_comparison_tests = load_file('next_gen_vs_2d_metaprogram_expression_tests.csv', local_dir=PROCESSED_DIR)
    mmp_volcano_df = mmp_comparison_tests[mmp_comparison_tests['OncotreeLineage'].str.contains("All Organoid Lineages")].copy()

    # add annotation columns to facilitate plotting
    mmp_volcano_df['-Mean Difference'] = -mmp_volcano_df['MeanDifference']
    mmp_volcano_df['MannWhitneyPValue'] = mmp_volcano_df['MannWhitneyPValue'].astype('float')
    mmp_volcano_df['p_transform'] = -np.log10(mmp_volcano_df['MannWhitneyPValue'])
    mmp_volcano_df['highlight'] = (
        ((mmp_volcano_df['-Mean Difference'] > mean_difference_threshold) | (mmp_volcano_df['MeanDifference'] > mean_difference_threshold)) &
        (mmp_volcano_df['MannWhitneyPValue'] < p_threshold)
    )

    plt.figure(figsize=(55*mm, 65*mm))
    plt.subplots_adjust(bottom=0.15, right=0.95, left=0.15)
    
    xlabel='-Mean Difference'
    ylabel='p_transform'
    
    sns.scatterplot(
        mmp_volcano_df,
        x=xlabel, y=ylabel, 
        hue='highlight', palette={True: highlight_color, False: 'tab:blue'},
        s=8,
        legend=False
    )
    plt.ylabel('-log10(p-value)', fontdict={'size': LABEL_SIZE})
    plt.xlabel('$\Delta$ Mean Metaprogram Score\n(2D - Organoid)', fontdict={'size': LABEL_SIZE})
    plt.title('2D vs. Organoid Metaprogram Scores', fontdict={'size': TITLE_SIZE})

    # add visual markers of thresholds
    plt.axhline(-np.log10(p_threshold), linestyle='dashed', color='tab:red', alpha=0.5, linewidth=0.5)
    plt.axvline(-mean_difference_threshold, linestyle='dashed', color='tab:red', alpha=0.5, linewidth=0.5)
    plt.axvline(mean_difference_threshold, linestyle='dashed', color='tab:red', alpha=0.5, linewidth=0.5)

    minx, maxx = plt.xlim()
    miny, maxy = plt.ylim()

    plt.gca().fill_between([minx, maxx], y1=[miny, miny], y2=[15, 15], facecolor='tab:gray', alpha=.1)
    plt.gca().fill_between([-0.5, 0.5], y1=[miny, miny], y2=[maxy, maxy], facecolor='tab:gray', alpha=.1)

    plt.xlim(minx, maxx)
    plt.ylim(miny, maxy)

    # annotate specific points with metaprogram names
    annots = mmp_volcano_df.query('highlight').sort_values('p_transform', ascending=False).set_index('Metaprogram').rename(
        {'Stress ': 'Stress', 'Cell Cycle - G2/M': 'G2/M', 'EMT-III ': 'EMT-III', 'Cell Cycle - G1/S': 'G1/S'}
    )

    manually_annotate(plt.gca(), "PDAC-classical", annots, 
                      xlabel, ylabel, (-0.1, -2), ha='left', va='center',
                      arrowprops=dict(arrowstyle='-', color='black', alpha=0, linewidth=0.25, shrinkA=0, shrinkB=1))
    manually_annotate(plt.gca(), "Secreted II", annots, 
                      xlabel, ylabel, (0, -1), ha='left', va='center',
                      arrowprops=dict(arrowstyle='-', color='black', alpha=0, linewidth=0.25, shrinkA=0, shrinkB=1))
    manually_annotate(plt.gca(), "Stress", annots, 
                      xlabel, ylabel, (0, -2), ha='left', va='center',
                      arrowprops=dict(arrowstyle='-', color='black', alpha=0, linewidth=0.25, shrinkA=0, shrinkB=1))
    manually_annotate(plt.gca(), "Cilia", annots, 
                      xlabel, ylabel, (0, -2), ha='right', va='center',
                      arrowprops=dict(arrowstyle='-', color='black', alpha=0, linewidth=0.25, shrinkA=0, shrinkB=0.5))
    manually_annotate(plt.gca(), "G2/M", annots, 
                      xlabel, ylabel, (-0.05, 2), ha='center', va='center',
                      arrowprops=dict(arrowstyle='-', color='black', alpha=0, linewidth=0.25, shrinkA=0, shrinkB=0))

    plt.savefig(os.path.join(FIGURE_DIR, 'Fig_4a_metaprogram_organoid_vs_adherent_volcano.pdf'))
    
# figure 4b

def plot_pdac_distribution_boxplots():
    """
    Plot 1D distributions of PDAC-classical expression across NextGen lineages between Depmap models and lineage-matched TCGA samples.
    """
    all_mmp_scores = load_file('celligner_geneset_scores.csv', local_dir=PROCESSED_DIR, index_col=0)
    omics_models_meta = load_file('model_metadata.csv', index_col=0).replace(lineage_replacement)
    all_celligner_samples = load_file('celligner_coordinates_w_hcmi.csv', index_col=0).replace(lineage_replacement)
    # lineages_to_show = ['Pancreas', 'Esophagus/Stomach', 'Bowel', 'Prostate', 'Breast', 'Ovary/Fallopian Tube']
    lineages_to_show = ['Pancreas', 'Esophagus/Stomach', 'Colorectal', 'Prostate', 'Breast', 'Ovary/Fallopian Tube', 'All Organoid Lineages']
    all_organoid_lineages = ['Esophagus/Stomach', 'Pancreas', 'Biliary Tract', 'CNS/Brain', 'Breast', 
                             'Colorectal', 'Bowel', 'Prostate', 'Ovary/Fallopian Tube', 'Ampulla of Vater', 'Uterus']

    # all expression scores per lineage, annotated by their respective lineages
    tirosh_distro_comparison_df = pd.concat([
        omics_models_meta[
            omics_models_meta['OncotreeLineage'].isin(lineages_to_show) & 
            omics_models_meta['IsNextGen'] & 
            omics_models_meta['GrowthPattern'].isin(['Dome'])
        ].assign(SampleSet = "DepMap Organoid")[['SampleSet', 'OncotreeLineage']],
        omics_models_meta[
            omics_models_meta['OncotreeLineage'].isin(lineages_to_show) & 
            omics_models_meta['IsAdherent2D']
        ].assign(SampleSet = "DepMap 2D")[['SampleSet', 'OncotreeLineage']],
        all_celligner_samples[
            all_celligner_samples['lineage'].isin(lineages_to_show) & 
            (all_celligner_samples['type'] == "TCGA+ tumor")
        ].assign(SampleSet = "TCGA+").rename({'lineage': 'OncotreeLineage'}, axis=1)[['SampleSet', 'OncotreeLineage']]
    ], axis=0)

    # all expression scores in any lineage, annotated by an aggregate label covering all lineages
    tirosh_overcounted_comparison_df = pd.concat([
        omics_models_meta[
            omics_models_meta['OncotreeLineage'].isin(all_organoid_lineages) & 
            omics_models_meta['IsNextGen'] & 
            omics_models_meta['GrowthPattern'].isin(['Dome'])
        ].assign(SampleSet = "DepMap Organoid")[['SampleSet', 'OncotreeLineage']],
        omics_models_meta[
            omics_models_meta['OncotreeLineage'].isin(all_organoid_lineages) & 
            omics_models_meta['IsAdherent2D']
        ].assign(SampleSet = "DepMap 2D")[['SampleSet', 'OncotreeLineage']],
        all_celligner_samples[
            all_celligner_samples['lineage'].isin(all_organoid_lineages) & 
            (all_celligner_samples['type'] == "TCGA+ tumor")
        ].assign(SampleSet = "TCGA+").rename({'lineage': 'OncotreeLineage'}, axis=1)[['SampleSet', 'OncotreeLineage']]
    ], axis=0)
    tirosh_overcounted_comparison_df['OncotreeLineage'] = 'All Organoid Lineages'

    # Run all hypothesis tests for the difference in expression score between models and tumors
    organoid_vs_tcga_tests = []
    adherent_vs_tcga_tests = []
    for lin in lineages_to_show:
        if lin != 'All Organoid Lineages':
            organoid_vs_tcga_tests.append(
                run_mann_whitney_u(
                    all_mmp_scores[['PDAC-classical']], 
                    tirosh_distro_comparison_df[
                        (tirosh_distro_comparison_df['SampleSet'] == 'DepMap Organoid') &
                        (tirosh_distro_comparison_df['OncotreeLineage'] == lin)
                    ].index.tolist(),
                    tirosh_distro_comparison_df[
                        (tirosh_distro_comparison_df['SampleSet'] == 'TCGA+') &
                        (tirosh_distro_comparison_df['OncotreeLineage'] == lin)
                    ].index.tolist(),
                    nan_policy='omit',
                    group2_name='TCGA+',
                ).assign(OncotreeLineage=lin)
            )
    
            adherent_vs_tcga_tests.append(
                run_mann_whitney_u(
                    all_mmp_scores[['PDAC-classical']], 
                    tirosh_distro_comparison_df[
                        (tirosh_distro_comparison_df['SampleSet'] == 'DepMap 2D') &
                        (tirosh_distro_comparison_df['OncotreeLineage'] == lin)
                    ].index.tolist(),
                    tirosh_distro_comparison_df[
                        (tirosh_distro_comparison_df['SampleSet'] == 'TCGA+') &
                        (tirosh_distro_comparison_df['OncotreeLineage'] == lin)
                    ].index.tolist(),
                    nan_policy='omit',
                    group2_name='TCGA+',
                ).assign(OncotreeLineage=lin)
            )
        else:
            organoid_vs_tcga_tests.append(
                run_mann_whitney_u(
                    all_mmp_scores[['PDAC-classical']], 
                    tirosh_distro_comparison_df[
                        (tirosh_distro_comparison_df['SampleSet'] == 'DepMap Organoid') &
                        (tirosh_distro_comparison_df['OncotreeLineage'].isin(all_organoid_lineages))
                    ].index.tolist(),
                    tirosh_distro_comparison_df[
                        (tirosh_distro_comparison_df['SampleSet'] == 'TCGA+') &
                        (tirosh_distro_comparison_df['OncotreeLineage'].isin(all_organoid_lineages))
                    ].index.tolist(),
                    nan_policy='omit',
                    group2_name='TCGA+',
                ).assign(OncotreeLineage=lin)
            )
    
            adherent_vs_tcga_tests.append(
                run_mann_whitney_u(
                    all_mmp_scores[['PDAC-classical']], 
                    tirosh_distro_comparison_df[
                        (tirosh_distro_comparison_df['SampleSet'] == 'DepMap 2D') &
                        (tirosh_distro_comparison_df['OncotreeLineage'].isin(all_organoid_lineages))
                    ].index.tolist(),
                    tirosh_distro_comparison_df[
                        (tirosh_distro_comparison_df['SampleSet'] == 'TCGA+') &
                        (tirosh_distro_comparison_df['OncotreeLineage'].isin(all_organoid_lineages))
                    ].index.tolist(),
                    nan_policy='omit',
                    group2_name='TCGA+',
                ).assign(OncotreeLineage=lin)
            )
    
    organoid_vs_tcga_tests = pd.concat(organoid_vs_tcga_tests, axis=0, ignore_index=True).set_index('OncotreeLineage')
    adherent_vs_tcga_tests = pd.concat(adherent_vs_tcga_tests, axis=0, ignore_index=True).set_index('OncotreeLineage')

    # join all data points into one table
    annotated_pdac_scores = all_mmp_scores[['PDAC-classical']].join(tirosh_distro_comparison_df, how='inner')
    overcounted_annotated_pdac_scores = all_mmp_scores[['PDAC-classical']].join(tirosh_overcounted_comparison_df, how='inner')
    all_annotated_pdac_scores = pd.concat([annotated_pdac_scores, overcounted_annotated_pdac_scores], axis=0)

    plt.figure(figsize=(75 * mm, 65 * mm))
    plt.subplots_adjust(left=0.25, right=0.95, bottom=0.15)

    # plot boxplots
    sns.boxplot(
        all_annotated_pdac_scores,
        x='PDAC-classical',
        y='OncotreeLineage',
        hue='SampleSet',
        palette=sample_set_palette,
        hue_order=['TCGA+', 'DepMap Organoid', 'DepMap 2D'],
        order=lineages_to_show,
        flierprops={'markersize': 2, 'markeredgewidth':0.25},
        whiskerprops={'linewidth': 0.5},
        boxprops={'linewidth': 0.5},
        medianprops={'linewidth': 0.5},
        capprops={'linewidth': 0.5}
    )

    xmin, xmax = plt.xlim()
    extend_xmax = xmax + (0.08 * (xmax - xmin))
    plt.xlim(xmin, extend_xmax)

    # plot significance bars
    for i, lin in enumerate(lineages_to_show):
        x_max_val = all_annotated_pdac_scores[all_annotated_pdac_scores['OncotreeLineage'] == lin]['PDAC-classical'].dropna().max()
        horizontal_annotate_significance_bar(
            plt.gca(), format_significance(organoid_vs_tcga_tests.loc[lin, 'FDR']), 
            x_max_val + 0.25, [i-1/3, i], text_offset=0.15, cap_length=-0.1 
        )
        horizontal_annotate_significance_bar(
            plt.gca(), format_significance(adherent_vs_tcga_tests.loc[lin, 'FDR']), 
            x_max_val + 0.75, [i-1/3, i+1/3], text_offset=0.15, cap_length=-0.1 
        )

    plt.xticks(plt.gca().get_xticks()[1:-1], labels=plt.gca().get_xticklabels()[1:-1], fontdict={'size': TICK_SIZE})
    plt.xlabel('PDAC-classical/EED Score', fontdict={'size': LABEL_SIZE})
    plt.yticks(plt.gca().get_yticks(), labels=lineages_to_show, rotation=0, verticalalignment='center', fontdict={'size': TICK_SIZE})
    plt.ylabel('')
    plt.title('PDAC-classical/EED Scores by Lineage')

    plt.legend(loc='lower left', bbox_to_anchor=(0.62, 0.39), ncol=1,
               handlelength=0.75, handletextpad=0.5, columnspacing=2, prop={'size': ANNOT_SIZE})
    plt.savefig(os.path.join(FIGURE_DIR, 'Fig_4b_pdac_distributions_dataset_comparison_by_lineage_box.pdf'))
    annotated_pdac_scores.rename_axis('SampleID').to_csv(os.path.join(PROCESSED_DIR, 'pdac_score_comparison_by_lineage.csv'))
    
    
def main():
    """
    Generate all plots related to metaprogram expression (gene sets from Gavish et al., 2023)
    """
    # figure 3a
    plot_metaprogram_bubble()
    
    # figure 4a
    plot_metaprogram_volcano(mean_difference_threshold=0.5, p_threshold=1e-30)
    
    # figure 4b
    plot_pdac_distribution_boxplots()

if __name__ == "__main__":
    main()    