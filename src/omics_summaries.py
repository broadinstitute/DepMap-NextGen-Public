import os
import pandas as pd
import numpy as np
from figure_utils import *
from constants import *
from data_utils import *

import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.patches import Patch, Rectangle
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.gridspec as gridspec

mpl.style.use(os.path.join(ASSETS_DIR, 'stylesheet.mplstyle'))

def load_data():
    """
    Load all data related to omics model summaries
    """
    omics_models_meta = load_file('model_metadata.csv', index_col=0)
    screen_metadata = load_file('screen_metadata.csv', index_col=0)
    cn = load_full_matrix('copy_number')
    damaging = load_full_matrix('damaging')
    hotspot = load_full_matrix('hotspot')
    omics_signatures = load_full_matrix('omics_signatures')
    alteration_frequency_df = load_file('selected_genomic_alteration_frequencies.csv', local_dir=PROCESSED_DIR)
    
    return omics_models_meta, screen_metadata, cn, damaging, hotspot, omics_signatures, alteration_frequency_df

# figure 1b

def make_omics_wheels(omics_models_meta):
    """
    Plot pie charts indicating which data modalities for NextGen models are covered by our dataset

    Args:
        omics_models_meta (pandas.DataFrame): A table of model metadata
    """
    next_gen_omics_metadata = omics_models_meta[omics_models_meta['IsNextGen']]
    
    fig, axs = plt.subplots(1, 3, figsize=(72*mm, 24*mm), gridspec_kw={'wspace': 0.3})
    plt.subplots_adjust(left=0.08, right=0.92)

    # Common parameters for wheel subplots
    r = 1.5
    wheel_width = 0.7
    title_y = 0.57
    n_y = title_y - 0.16

    # WGS plot
    wgs_wheel = next_gen_omics_metadata[next_gen_omics_metadata['HasWGS']].value_counts('OncotreeLineage')
    axs[0].pie(
        wgs_wheel, radius=r, wedgeprops=dict(width=wheel_width), counterclock=False, startangle=90, colors=[lineage_cmap[lin] for lin in wgs_wheel.index.tolist()]
    )
    axs[0].text(0.5, title_y, 'WGS', transform=axs[0].transAxes, horizontalalignment='center', verticalalignment='center', 
                fontdict={'size': TITLE_SIZE})
    axs[0].text(0.5, n_y, f'(n = {wgs_wheel.sum()})', transform=axs[0].transAxes, horizontalalignment='center', verticalalignment='center',
                fontdict={'size': LABEL_SIZE})

    # RNA plot
    rna_wheel = next_gen_omics_metadata[next_gen_omics_metadata['HasExprData']].value_counts('OncotreeLineage')
    axs[1].pie(
        rna_wheel, radius=r, wedgeprops=dict(width=wheel_width), counterclock=False, startangle=90, colors=[lineage_cmap[lin] for lin in rna_wheel.index.tolist()]
    )
    axs[1].text(0.5, title_y, 'RNAseq', transform=axs[1].transAxes, horizontalalignment='center', verticalalignment='center', 
                fontdict={'size': TITLE_SIZE})
    axs[1].text(0.5, n_y, f'(n = {rna_wheel.sum()})', transform=axs[1].transAxes, horizontalalignment='center', verticalalignment='center',
                fontdict={'size': LABEL_SIZE})

    # CRISPR plot
    crispr_wheel = next_gen_omics_metadata[next_gen_omics_metadata['CRISPRScreenType'].isin(['2DO', '3DO', '3DN', '2DN'])].value_counts('OncotreeLineage')
    axs[2].pie(
        crispr_wheel, radius=r, wedgeprops=dict(width=wheel_width), counterclock=False, startangle=90, colors=[lineage_cmap[lin] for lin in crispr_wheel.index.tolist()]
    )
    axs[2].text(0.5, title_y, 'CRISPR', transform=axs[2].transAxes, horizontalalignment='center', verticalalignment='center', 
                fontdict={'size': TITLE_SIZE})
    axs[2].text(0.5, n_y, f'(n = {crispr_wheel.sum()})', transform=axs[2].transAxes, horizontalalignment='center', verticalalignment='center',
                fontdict={'size': LABEL_SIZE})

    plt.savefig(os.path.join(FIGURE_DIR, f'Fig_1b_omics_wheels.pdf'))


# figure 1c

def prepare_paneled_heatmap_axes(figsize=(72 * mm, 30 * mm), 
                                 column_ratios=[1, 1], column_widthspace=0.5,
                                 row_ratios = [[1, 1], [1, 1, 1, 1]], row_heightspace=[0.5, 0.5],
                                 cbar_ratios = [[9, 1], [3, 1]], cbar_heightspace=[0.05, 0.1]):
    """
    Create figure axes for plotting frequent gene alterations in respective lineages
    """
    fig = plt.figure(figsize=figsize)
    # two columns
    outer = gridspec.GridSpec(1, 2, width_ratios=column_ratios, wspace=column_widthspace)
    # left column, two rows
    gsl = gridspec.GridSpecFromSubplotSpec(subplot_spec=outer[0], nrows=2, ncols=1, hspace=row_heightspace[0], height_ratios=row_ratios[0])
    # left: top and bottom
    gslt = gridspec.GridSpecFromSubplotSpec(subplot_spec=gsl[0], nrows=2, ncols=1, hspace=cbar_heightspace[0], height_ratios=cbar_ratios[0])
    gslb = gridspec.GridSpecFromSubplotSpec(subplot_spec=gsl[1], nrows=2, ncols=1, hspace=cbar_heightspace[0], height_ratios=cbar_ratios[0])
    # right column, four rows
    gsr = gridspec.GridSpecFromSubplotSpec(subplot_spec=outer[1], nrows=4, ncols=1, hspace=row_heightspace[1], height_ratios=row_ratios[1])
    # right: top, 2nd, 3rd, and bottom
    gsrt = gridspec.GridSpecFromSubplotSpec(subplot_spec=gsr[0], nrows=2, ncols=1, hspace=cbar_heightspace[1], height_ratios=cbar_ratios[1])
    gsr2 = gridspec.GridSpecFromSubplotSpec(subplot_spec=gsr[1], nrows=2, ncols=1, hspace=cbar_heightspace[1], height_ratios=cbar_ratios[1])
    gsr3 = gridspec.GridSpecFromSubplotSpec(subplot_spec=gsr[2], nrows=2, ncols=1, hspace=cbar_heightspace[1], height_ratios=cbar_ratios[1])
    gsrb = gridspec.GridSpecFromSubplotSpec(subplot_spec=gsr[3], nrows=2, ncols=1, hspace=cbar_heightspace[1], height_ratios=cbar_ratios[1])

    ax_map = dict()
    # left
    ax_map[0] = dict()
    # left top
    ax_map[0][0] = dict()
    ax_map[0][0][0] = fig.add_subplot(gslt[0])
    ax_map[0][0][1] = fig.add_subplot(gslt[1])
    # left bottom
    ax_map[0][1] = dict()
    ax_map[0][1][0] = fig.add_subplot(gslb[0])
    ax_map[0][1][1] = fig.add_subplot(gslb[1])
    # right
    ax_map[1] = dict()
    # right top
    ax_map[1][0] = dict()
    ax_map[1][0][0] = fig.add_subplot(gsrt[0])
    ax_map[1][0][1] = fig.add_subplot(gsrt[1])
    # right 2nd
    ax_map[1][1] = dict()
    ax_map[1][1][0] = fig.add_subplot(gsr2[0])
    ax_map[1][1][1] = fig.add_subplot(gsr2[1])
    # right 3rd
    ax_map[1][2] = dict()
    ax_map[1][2][0] = fig.add_subplot(gsr3[0])
    ax_map[1][2][1] = fig.add_subplot(gsr3[1])
    # right bottom
    ax_map[1][3] = dict()
    ax_map[1][3][0] = fig.add_subplot(gsrb[0])
    ax_map[1][3][1] = fig.add_subplot(gsrb[1])
    
    return fig, ax_map
    
def make_alteration_heatmap_panels(alteration_heatmap_df):
    """
    Plot alteration frequencies per lineage and dataset for commonly mutated genes

    Args:
        alteration_heatmap_df (pandas.DataFrame): A table of alteration frequencies per lineage and dataset
    """
    fig, ax_map = prepare_paneled_heatmap_axes(figsize=(72 * mm, 48 * mm), column_widthspace=2, row_heightspace=[0.4, 1])

    common_freq_values = [0.4, 0.6, 0.8, 1]
    common_freq_labels = ['40', '60', '80', '      100 (%)']
    rare_freq_values = [0, 0.05, 0.1, 0.15]
    rare_freq_labels = ['0', '5', '10', '      15 (%)']


    # upper left
    common_colorectal_matrix = alteration_heatmap_df[alteration_heatmap_df['OncotreeLineage'] == "Bowel"].pivot(
        index='FeatureLabel', columns='Dataset', values='Frequency'
    ).reindex(
        index=['TP53 damaging', 'APC damaging', '18q loss'], columns=['Tumor', 'DepMap NextGen', 'DepMap Traditional']
    )
    lcm = LinearSegmentedColormap.from_list('Bowel', ['white', lineage_cmap['Bowel']])
    lnorm = Normalize(vmin=common_freq_values[0], vmax=common_freq_values[-1])
    sns.heatmap(
        common_colorectal_matrix, ax=ax_map[0][0][0], cbar=False, cmap=lcm, norm=lnorm, linewidths=0.25, linecolor='black'
    )
    ax_map[0][0][0].set_xlabel('')
    ax_map[0][0][0].set_ylabel('')
    ax_map[0][0][0].xaxis.set_ticks_position('top')
    ax_map[0][0][0].set_xticklabels(labels=common_colorectal_matrix.columns.tolist(), rotation=90)
    ax_map[0][0][0].text(-0.25, -0.25, 'Colorectal', fontdict={'fontsize': ANNOT_SIZE, 'fontweight': 'bold'}, ha='right')

    sm = plt.cm.ScalarMappable(cmap=lcm, norm=lnorm)
    sm.set_array([])
    cbar = plt.colorbar(sm, cax=ax_map[0][0][1], orientation='horizontal')
    cbar.set_ticks(common_freq_values, labels=common_freq_labels, fontsize=TICK_SIZE)
    cbar.ax.tick_params(width=0.25)
    cbar.outline.set_linewidth(0.25)


    # lower left
    common_pancreas_matrix = alteration_heatmap_df[alteration_heatmap_df['OncotreeLineage'] == "Pancreas"].pivot(
        index='FeatureLabel', columns='Dataset', values='Frequency'
    ).reindex(
        index=['TP53 damaging', 'KRAS hotspot', '9p loss'], columns=['Tumor', 'DepMap NextGen', 'DepMap Traditional']
    )
    lcm = LinearSegmentedColormap.from_list('Pancreas', ['white', lineage_cmap['Pancreas']])
    lnorm = Normalize(vmin=common_freq_values[0], vmax=common_freq_values[-1])
    sns.heatmap(
        common_pancreas_matrix, ax=ax_map[0][1][0], cbar=False, cmap=lcm, norm=lnorm, linewidths=0.25, linecolor='black'
    )
    ax_map[0][1][0].set_xlabel('')
    ax_map[0][1][0].set_ylabel('')
    ax_map[0][1][0].set_xticks([])
    ax_map[0][1][0].text(-0.25, -0.25, 'Pancreas', fontdict={'fontsize': ANNOT_SIZE, 'fontweight': 'bold'}, ha='right')

    sm = plt.cm.ScalarMappable(cmap=lcm, norm=lnorm)
    sm.set_array([])
    cbar = plt.colorbar(sm, cax=ax_map[0][1][1], orientation='horizontal')
    cbar.set_ticks(common_freq_values, labels=common_freq_labels, fontsize=TICK_SIZE)
    cbar.ax.tick_params(width=0.25)
    cbar.outline.set_linewidth(0.25)


    # upper right
    common_cns_matrix = alteration_heatmap_df[alteration_heatmap_df['OncotreeLineage'] == "CNS/Brain"].pivot(
        index='FeatureLabel', columns='Dataset', values='Frequency'
    ).reindex(
        index=['TERT promoter hotspot'], columns=['Tumor', 'DepMap NextGen', 'DepMap Traditional']
    )
    lcm = LinearSegmentedColormap.from_list('CNS/Brain', ['white', lineage_cmap['CNS/Brain']])
    lnorm = Normalize(vmin=common_freq_values[0], vmax=common_freq_values[-1])
    sns.heatmap(
        common_cns_matrix, ax=ax_map[1][0][0], cbar=False, cmap=lcm, norm=lnorm, linewidths=0.25, linecolor='black'
    )
    ax_map[1][0][0].set_xlabel('')
    ax_map[1][0][0].set_ylabel('')
    ax_map[1][0][0].xaxis.set_ticks_position('top')
    ax_map[1][0][0].set_xticklabels(labels=common_cns_matrix.columns.tolist(), rotation=90)
    ax_map[1][0][0].set_yticklabels(labels=['TERT promoter hotspot'], rotation=0)
    ax_map[1][0][0].text(-0.25, -0.25, 'CNS', fontdict={'fontsize': ANNOT_SIZE, 'fontweight': 'bold'}, ha='right')

    sm = plt.cm.ScalarMappable(cmap=lcm, norm=lnorm)
    sm.set_array([])
    cbar = plt.colorbar(sm, cax=ax_map[1][0][1], orientation='horizontal')
    cbar.set_ticks(common_freq_values, labels=common_freq_labels, fontsize=TICK_SIZE)
    cbar.ax.tick_params(width=0.25)
    cbar.outline.set_linewidth(0.25)


    # 2nd right
    rare_colorectal_matrix = alteration_heatmap_df[alteration_heatmap_df['OncotreeLineage'] == "Bowel"].pivot(
        index='FeatureLabel', columns='Dataset', values='Frequency'
    ).reindex(
        index=['NRAS hotspot'], columns=['Tumor', 'DepMap NextGen', 'DepMap Traditional']
    )
    lcm = LinearSegmentedColormap.from_list('Bowel', ['white', lineage_cmap['Bowel']])
    lnorm = Normalize(vmin=rare_freq_values[0], vmax=rare_freq_values[-1])
    sns.heatmap(
        rare_colorectal_matrix, ax=ax_map[1][1][0], cbar=False, cmap=lcm, norm=lnorm, linewidths=0.25, linecolor='black'
    )
    ax_map[1][1][0].set_xlabel('')
    ax_map[1][1][0].set_ylabel('')
    ax_map[1][1][0].set_xticks([])
    ax_map[1][1][0].set_yticklabels(labels=['NRAS hotspot'], rotation=0)
    ax_map[1][1][0].text(-0.25, -0.25, 'Colorectal', fontdict={'fontsize': ANNOT_SIZE, 'fontweight': 'bold'}, ha='right')

    sm = plt.cm.ScalarMappable(cmap=lcm, norm=lnorm)
    sm.set_array([])
    cbar = plt.colorbar(sm, cax=ax_map[1][1][1], orientation='horizontal')
    cbar.set_ticks(rare_freq_values, labels=rare_freq_labels, fontsize=TICK_SIZE)
    cbar.ax.tick_params(width=0.25)
    cbar.outline.set_linewidth(0.25)


    # 3rd right
    rare_breast_matrix = alteration_heatmap_df[alteration_heatmap_df['OncotreeLineage'] == "Breast"].pivot(
        index='FeatureLabel', columns='Dataset', values='Frequency'
    ).reindex(
        index=['ESR1 hotspot'], columns=['Tumor', 'DepMap NextGen', 'DepMap Traditional']
    )
    lcm = LinearSegmentedColormap.from_list('Breast', ['white', lineage_cmap['Breast']])
    lnorm = Normalize(vmin=rare_freq_values[0], vmax=rare_freq_values[-1])
    sns.heatmap(
        rare_breast_matrix, ax=ax_map[1][2][0], cbar=False, cmap=lcm, norm=lnorm, linewidths=0.25, linecolor='black'
    )
    ax_map[1][2][0].set_xlabel('')
    ax_map[1][2][0].set_ylabel('')
    ax_map[1][2][0].set_xticks([])
    ax_map[1][2][0].set_yticklabels(labels=['ESR1 hotspot'], rotation=0)
    ax_map[1][2][0].text(-0.25, -0.25, 'Breast', fontdict={'fontsize': ANNOT_SIZE, 'fontweight': 'bold'}, ha='right')

    sm = plt.cm.ScalarMappable(cmap=lcm, norm=lnorm)
    sm.set_array([])
    cbar = plt.colorbar(sm, cax=ax_map[1][2][1], orientation='horizontal')
    cbar.set_ticks(rare_freq_values, labels=rare_freq_labels, fontsize=TICK_SIZE)
    cbar.ax.tick_params(width=0.25)
    cbar.outline.set_linewidth(0.25)


    # bottom right
    rare_prostate_matrix = alteration_heatmap_df[alteration_heatmap_df['OncotreeLineage'] == "Prostate"].pivot(
        index='FeatureLabel', columns='Dataset', values='Frequency'
    ).reindex(
        index=['SPOP damaging'], columns=['Tumor', 'DepMap NextGen', 'DepMap Traditional']
    )
    lcm = LinearSegmentedColormap.from_list('Prostate', ['white', lineage_cmap['Prostate']])
    lnorm = Normalize(vmin=rare_freq_values[0], vmax=rare_freq_values[-1])
    sns.heatmap(
        rare_prostate_matrix, ax=ax_map[1][3][0], cbar=False, cmap=lcm, norm=lnorm, linewidths=0.25, linecolor='black'
    )
    ax_map[1][3][0].set_xlabel('')
    ax_map[1][3][0].set_ylabel('')
    ax_map[1][3][0].set_xticks([])
    ax_map[1][3][0].set_yticklabels(labels=['SPOP damaging'], rotation=0)
    ax_map[1][3][0].text(-0.25, -0.25, 'Prostate', fontdict={'fontsize': ANNOT_SIZE, 'fontweight': 'bold'}, ha='right')

    sm = plt.cm.ScalarMappable(cmap=lcm, norm=lnorm)
    sm.set_array([])
    cbar = plt.colorbar(sm, cax=ax_map[1][3][1], orientation='horizontal')
    cbar.set_ticks(rare_freq_values, labels=rare_freq_labels, fontsize=TICK_SIZE)
    cbar.ax.tick_params(width=0.25)
    cbar.outline.set_linewidth(0.25)

    # add publication markings
    ax_map[0][0][0].add_patch(Rectangle((0, 2), 1, 1, fill=False, edgecolor='red', lw=0.5, clip_on=False))
    ax_map[0][1][0].add_patch(Rectangle((0, 2), 1, 1, fill=False, edgecolor='red', lw=0.5, clip_on=False))
    ax_map[1][0][0].add_patch(Rectangle((0, 0), 1, 1, fill=False, edgecolor='red', lw=0.5, clip_on=False))
    ax_map[1][1][0].add_patch(Rectangle((0, 0), 1, 1, fill=False, edgecolor='red', lw=0.5, clip_on=False))
    ax_map[1][2][0].add_patch(Rectangle((0, 0), 1, 1, fill=False, edgecolor='red', lw=0.5, clip_on=False))
    ax_map[1][3][0].add_patch(Rectangle((0, 0), 1, 1, fill=False, edgecolor='red', lw=0.5, clip_on=False))

    plt.subplots_adjust(right=0.92, left=0.24, top=0.65, bottom=0.1)
    plt.savefig(os.path.join(FIGURE_DIR, 'Fig_1c_prominent_alteration_heatmaps.pdf'))
    

# ED 1b

def make_omics_summary(mutations, omics_models_meta, screen_metadata, cn, hotspot, damaging, omics_signatures):
    """
    Plot an annotated heatmap across all NextGen models showing broad phenotypic patterns

    Args:
        mutations (list[str]):  A list of genes to show mutations for
        omics_models_meta (pandas.DataFrame): A table of model metadata
        screen_metadata (pandas.DataFrame): A table of screen metadata
        cn (pandas.DataFrame): A matrix (model x gene) of relative copy number
        hotspot (pandas.DataFrame): A boolean matrix (model x gene) indicating the presence of hotspot mutations
        damaging (pandas.DataFrame): A boolean matrix (model x gene) indicating the presence of damaging mutations
        omics_signatures (pandas.DataFrame): A matrix (model x phenotype) of phenotype values
    """
    # Define models to use
    omics_summary_models = omics_models_meta[
        omics_models_meta['IsNextGen'] &
        ~omics_models_meta.index.isin(screen_metadata[screen_metadata['ScreenType'] == "2DO"]['ModelID'].tolist())
    ]

    # Define genes to show in mutation panel
    mutation_list = [search_gene(g) for g in mutations]

    # Define the mutation matrix by merging the hotspot and damaging matrices
    h, d = hotspot.align(damaging, join='outer')
    hotspot_or_damaging = ((h > 0) | (d > 0)).reindex(columns=mutation_list)

    # Define the cn matrix, sorting genes by cytoband
    cn_genes_sorted_metadata = hgnc_annotations[
        hgnc_annotations['cds_gene_id'].isin(cn.columns.tolist())
    ].sort_values('location_sortable', ignore_index=True).dropna(subset=['location_sortable'])
    cn_genes_sorted_metadata = cn_genes_sorted_metadata[
        cn_genes_sorted_metadata['cds_gene_id'].isin(cn.var().loc[lambda x: x > 0.15].index.tolist())
    ].reset_index(drop=True)
    cn_genes_sorted = cn_genes_sorted_metadata['cds_gene_id'].tolist()
    cn_genome_sorted = cn.reindex(columns=cn_genes_sorted).dropna(how='all', axis=1)
    cn_transform_genome_sorted = np.log2(cn_genome_sorted+1).add_prefix('CN-')

    center_9p = (cn_genes_sorted_metadata[cn_genes_sorted_metadata['location_sortable'].str.startswith('09p')].index.min() + 
                 cn_genes_sorted_metadata[cn_genes_sorted_metadata['location_sortable'].str.startswith('09p')].index.max()) / 2
    center_18q = (cn_genes_sorted_metadata[cn_genes_sorted_metadata['location_sortable'].str.startswith('18q')].index.min() + 
                  cn_genes_sorted_metadata[cn_genes_sorted_metadata['location_sortable'].str.startswith('18q')].index.max()) / 2

    # Define the screened boolean series
    is_screened = pd.DataFrame({
        'Screened': ~omics_models_meta['CRISPRScreenType'].isna() & (omics_models_meta['CRISPRScreenType'] != '2DO'),
    }, index=omics_models_meta.index.tolist()).fillna(False)

    # Combine matrices to figure out the row order
    order_matrix = pd.concat([
        omics_models_meta[['OncotreeLineage']],
        hotspot_or_damaging.reindex(columns=mutation_list),
        is_screened
    ], axis=1)

    # Define a categorical to set a custom order for lineages
    order_matrix['OncotreeLineage'] = pd.Categorical(
        order_matrix['OncotreeLineage'], 
        categories=[
            'CNS/Brain', 'Esophagus/Stomach', 'Pancreas', 'Biliary Tract', 'Bowel', 'Breast', 'Ovary/Fallopian Tube', 'Prostate', 'Ampulla of Vater', 'Uterus'
        ],
        ordered=True
    )

    # Sort the matrix along combined features
    order_matrix = order_matrix.reindex(index=omics_summary_models.index.tolist()).sort_values(
        ['OncotreeLineage'] + mutation_list + ['Screened'], 
        ascending = [True] + [False for _ in mutation_list] + [False]
    )
    row_order = order_matrix.index.tolist()

    # Identify the lineage midpoints for y-axis labels
    row_label_midpoints = (
        (omics_models_meta.reindex(index=row_order).loc[:, ['OncotreeLineage']].reset_index().reset_index().groupby(
            'OncotreeLineage').idxmin(numeric_only=True) +
         omics_models_meta.reindex(index=row_order).loc[: , ['OncotreeLineage']].reset_index().reset_index().groupby(
             'OncotreeLineage').idxmax(numeric_only=True)) / 2
    )['index'].to_dict()

    # Produce the omics summary plot
    fig, ax_map = prepare_gridspec(
        pd.DataFrame([[True, True, True, True, True]]),
        figsize = (72*mm, 90*mm),
        width_ratios=[0.5, len(mutation_list) * 2 / 3, 3, 2, 0.25], wspace=0.1,
        inner_gs_dict={
            (0, 3): {'nrows':1, 'ncols':3, 'wspace': 0},
            (0, 4): {'nrows':7, 'ncols':1, 'hspace': 0, 'height_ratios': [1, 4, 1, 4, 1, 4, 1]}
        }, 
    )
    plt.subplots_adjust(left=0.28, top=0.9, bottom=0.2, right=0.9)

    # Pad colorbars by using dummy axes that get removed
    ax_map[0][4][(0, 0)].set_visible(False)
    ax_map[0][4][(2, 0)].set_visible(False)
    ax_map[0][4][(4, 0)].set_visible(False)
    ax_map[0][4][(6, 0)].set_visible(False)

    # Lineage
    lin_heatmap = omics_models_meta.reindex(index=row_order).loc[: , ['OncotreeLineage']]
    lin_heatmap['OncotreeLineage'], lin_heatmap_enum, lin_heatmap_cmap = convert_categoricals_by_cmap(lin_heatmap['OncotreeLineage'], lineage_cmap)
    sns.heatmap(
        lin_heatmap, cmap=lin_heatmap_cmap, cbar=False,
        yticklabels=False, xticklabels=True, 
        ax=ax_map[0][0]
    )
    ax_map[0][0].set_ylabel('')
    ax_map[0][0].set_xticks([0.5], ['Lineage'], rotation=60, horizontalalignment='right', fontsize=TICK_SIZE)
    offset_x_labels(fig, ax_map[0][0], 5, 0)

    for k in row_label_midpoints:
        if k == 'Uterus':
            y_off = 5
        elif k == 'Prostate':
            y_off = -2
        elif k == 'Other':
            y_off = 8
        else:
            y_off = 0

        ax_map[0][0].annotate(k, xy=(0, row_label_midpoints[k] + 0.5), xycoords=ax_map[0][0].transData, 
                              xytext=(-1, row_label_midpoints[k] + 0.5 + y_off), arrowprops={'arrowstyle': '-', 'color': 'black', 'shrinkA': 0, 'shrinkB': 0, 'linewidth': 0.5}, 
                              annotation_clip=False, horizontalalignment='right', verticalalignment='center', fontsize=TICK_SIZE)

    # Mutations
    mutation_heatmap = hotspot_or_damaging.reindex(index=row_order).loc[:, mutation_list].replace({True: 1, False: 0})
    mutation_heatmap = mutation_heatmap.rename({g:g.split(' ')[0] for g in mutation_heatmap.columns}, axis=1)
    mutation_heatmap_cmap = ListedColormap({0: 'white', 1: 'black'}.values())
    sns.heatmap(
        mutation_heatmap, cmap=mutation_heatmap_cmap, ax=ax_map[0][1], yticklabels=False, xticklabels=True, cbar=False
    )
    ax_map[0][1].set_ylabel('')
    ax_map[0][1].set_facecolor('lightgray')
    ax_map[0][1].set_xticks([0.5 + i for i in range(len(mutation_list))], 
                            [g.split(' ')[0] for g in mutation_list], 
                            rotation=60, horizontalalignment='right', fontsize=TICK_SIZE)
    offset_x_labels(fig, ax_map[0][1], 5, 0)
    ax_map[0][1].set_title('Mutations', fontdict={'fontsize': TITLE_SIZE})

    # CN
    cn_heatmap = cn_transform_genome_sorted.reindex(index=row_order)
    sns.heatmap(
        cn_heatmap, cmap='RdBu_r', ax=ax_map[0][2], yticklabels=False, xticklabels=False, vmin=0, vmax=2, cbar=False
    )
    ax_map[0][2].set_ylabel('')
    ax_map[0][2].set_xticks([center_9p, center_18q], 
                            ['9p', '18q'], 
                            rotation=60, horizontalalignment='right', fontsize=TICK_SIZE)
    offset_x_labels(fig, ax_map[0][2], 5, 0)
    ax_map[0][2].set_facecolor('lightgray')
    ax_map[0][2].set_title('CN', fontdict={'fontsize': TITLE_SIZE})

    # MSI Score, Ploidy, Screened
    msi_heatmap = omics_signatures.reindex(index=row_order).loc[:, ['MSIScore']]
    sns.heatmap(
        msi_heatmap, cmap='Greens', ax=ax_map[0][3][(0, 0)], yticklabels=False, xticklabels=False, vmin=0, vmax=10, cbar=False
    )
    ax_map[0][3][(0, 0)].set_ylabel('')
    ax_map[0][3][(0, 0)].set_facecolor('lightgray')
    ax_map[0][3][(0, 0)].set_xticks([0.5], ['MSIScore'], rotation=60, horizontalalignment='right', fontsize=TICK_SIZE)
    offset_x_labels(fig, ax_map[0][3][(0, 0)], 5, 0)

    ploidy_heatmap = omics_signatures.reindex(index=row_order).loc[:, ['Ploidy']]
    sns.heatmap(
        ploidy_heatmap, cmap='Purples', ax=ax_map[0][3][(0, 1)], yticklabels=False, xticklabels=False, vmin=0, vmax=4, cbar=False
    )
    ax_map[0][3][(0, 1)].set_ylabel('')
    ax_map[0][3][(0, 1)].set_facecolor('lightgray')
    ax_map[0][3][(0, 1)].set_xticks([0.5], ['Ploidy'], rotation=60, horizontalalignment='right', fontsize=TICK_SIZE)
    offset_x_labels(fig, ax_map[0][3][(0, 1)], 5, 0)

    screened_heatmap = is_screened.reindex(index=row_order).fillna(False).replace({True: 1, False: 0})
    screened_heatmap_cmap = ListedColormap({0: 'white', 1: 'black'}.values())
    sns.heatmap(
        screened_heatmap, cmap=screened_heatmap_cmap, ax=ax_map[0][3][(0, 2)], yticklabels=False, xticklabels=False, cbar=False
    )
    ax_map[0][3][(0, 2)].set_ylabel('')
    ax_map[0][3][(0, 2)].set_xticks([0.5], ['Screened'], rotation=60, horizontalalignment='right', fontsize=TICK_SIZE)
    offset_x_labels(fig, ax_map[0][3][(0, 2)], 5, 0)

    # Colorbars
    sm = plt.cm.ScalarMappable(cmap='RdBu_r', norm=Normalize(vmin=0, vmax=2))
    sm.set_array([])
    cbar = plt.colorbar(sm, cax=ax_map[0][4][(1, 0)], orientation='vertical')
    cbar.set_ticks([0, 1, 2], labels=[0, 1, 2], fontsize=TICK_SIZE)
    cbar.set_label('log2(Relative CN + 1)', rotation=90, labelpad=3, fontsize=TICK_SIZE)
    cbar.ax.yaxis.set_label_position('right')
    cbar.outline.set_linewidth(0.25)

    sm = plt.cm.ScalarMappable(cmap='Greens', norm=Normalize(vmin=0, vmax=10))
    sm.set_array([])
    cbar = plt.colorbar(sm, cax=ax_map[0][4][(3, 0)], orientation='vertical')
    cbar.set_ticks([2, 5, 8], labels=[2, 5, 8], fontsize=TICK_SIZE)
    cbar.set_label('MSI Score', rotation=90, labelpad=3, fontsize=TICK_SIZE)
    cbar.ax.yaxis.set_label_position('right')
    cbar.outline.set_linewidth(0.25)

    sm = plt.cm.ScalarMappable(cmap='Purples', norm=Normalize(vmin=0, vmax=4))
    sm.set_array([])
    cbar = plt.colorbar(sm, cax=ax_map[0][4][(5, 0)], orientation='vertical')
    cbar.set_ticks([0, 2, 4], labels=[0, 2, 4], fontsize=TICK_SIZE)
    cbar.set_label('Ploidy', rotation=90, labelpad=3, fontsize=TICK_SIZE)
    cbar.ax.yaxis.set_label_position('right')
    cbar.outline.set_linewidth(0.25)

    plt.savefig(os.path.join(FIGURE_DIR, 'Fig_ED1b_omics_summary.pdf'))
    
    
def main():
    """
    Generate all plots related to broad/survey patterns within our dataset.
    """
    omics_models_meta, screen_metadata, cn, damaging, hotspot, omics_signatures, alteration_frequency_df = load_data()
    # figure 1b
    make_omics_wheels(omics_models_meta.replace({'Lung': 'Other', 'Head and Neck': 'Other'}))
    
    # figure 1c
    make_alteration_heatmap_panels(alteration_frequency_df)
    
    # figure ED1b
    mutation_list = ['TP53', 'KRAS', 'CDKN2A', 'PIK3CA', 'BRAF', 'CTNNB1', 'APC', 'RB1', 'TERT', 'PTEN']
    make_omics_summary(mutation_list, omics_models_meta.replace({'Lung': 'Other', 'Head and Neck': 'Other'}), screen_metadata, cn, damaging, hotspot, omics_signatures)

if __name__ == "__main__":
    main()