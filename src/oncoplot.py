import os
import pandas as pd
import numpy as np

from figure_utils import *
from constants import *
from data_utils import *

import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.patches import Patch, Rectangle
from matplotlib.lines import Line2D
from matplotlib.legend_handler import HandlerPatch
import matplotlib.gridspec as gridspec

mpl.style.use(os.path.join(ASSETS_DIR, 'stylesheet.mplstyle'))

DAMAGING_BLACKLIST = ['SOS1', 'IDH1', 'ERBB2', 'NFE2L2', 'ESR1', 'BRAF', 'ERBB3']

# custom legend elements

from assets.custom_legend_elements import ForwardSlashRectangle, BackSlashRectangle, XRectangle, DotRectangle, BackSlashHandler, ForwardSlashHandler, XHandler, DotRectangleHandler, LargeBackSlashHandler, LargeForwardSlashHandler, LargeXHandler, LargeDotRectangleHandler, custom_handler_map, large_custom_handler_map

# load data

def load_data():
    omics_models_meta = load_file('model_metadata.csv', index_col=0)
    screen_metadata = load_file('screen_metadata.csv', index_col=0)
    geneset_table = load_file('geneset_table_updated.csv')

    oncogenes = geneset_table[geneset_table['Geneset'] == 'Oncogenes']['cds_gene_id'].tolist()
    tumor_suppressors = geneset_table[geneset_table['Geneset'] == 'Tumor suppressors']['cds_gene_id'].tolist()

    cn = load_full_matrix('copy_number')
    damaging = load_full_matrix('damaging')
    hotspot = load_full_matrix('hotspot')

    # mask damaging matrix for oncogenes
    damaging = damaging.drop([search_gene(g) for g in DAMAGING_BLACKLIST], axis=1, errors='ignore')
    damaging.loc[:, damaging.columns.isin(oncogenes)] = 0
    
    mutation_matrix = ((hotspot > 0) | (damaging > 0))
    cn_transform = discretize_matrix_by_interval(
        cn,
        interval_dict={'Amplification': pd.Interval(left=3, right=np.inf, closed='neither'),
                       'Deletion': pd.Interval(left=-np.inf, right=0.25, closed='neither'),
                       'Neutral Copy': pd.Interval(left=0.25, right=3, closed='both')}
    )

    # mask LoF oncogene and GoF tumor suppressor alterations
    cn_transform.loc[
        :, cn_transform.columns.isin(oncogenes)
    ] = cn_transform.loc[:, cn_transform.columns.isin(oncogenes)].replace({'Deletion': 'Neutral Copy'})
    cn_transform.loc[
        :, cn_transform.columns.isin(tumor_suppressors)
    ] = cn_transform.loc[:, cn_transform.columns.isin(tumor_suppressors)].replace({'Amplification': 'Neutral Copy'})
    
    return omics_models_meta, screen_metadata, damaging, hotspot, mutation_matrix, cn_transform

# helper functions

def create_cmap_from_series(series, palette='pastel', as_hex=True):
    lut = dict(zip(series.unique().tolist(), sns.color_palette(palette=palette, n_colors=len(set(series)))))
    if as_hex:
        lut = {k:to_hex(v) for k,v in lut.items()}
    return lut

def discretize_matrix_by_interval(value_matrix, interval_dict):
    discretized_matrix = value_matrix.copy()
    discretized_matrix.loc[:, :] = np.nan
    replacements = dict()
    i = 0
    for lb, interval in interval_dict.items():
        left_boolean = value_matrix >= interval.left if interval.closed_left else value_matrix > interval.left
        right_boolean = value_matrix <= interval.right if interval.closed_left else value_matrix < interval.right
        in_interval = (left_boolean & right_boolean).astype(bool)
        discretized_matrix[in_interval] = i
        replacements[i] = lb
        i += 1
    return discretized_matrix.replace(replacements)

def generate_indicator_sequence_map(value_encoding):
    df_length = int(np.power(2, len(value_encoding)))
    boolean_matrix = pd.DataFrame(columns=list(value_encoding.keys()), index=range(df_length))
    current_indicator_sequence = pd.Series(False, index=list(value_encoding.keys()))
    counters = {lb: 0 for lb in list(value_encoding.keys())}
    for i in range(df_length):
        boolean_matrix.loc[i, :] = current_indicator_sequence
        counters = {lb: counters[lb] + 1 for lb in counters}
        for lb, cnt in counters.items():
            if cnt == value_encoding[lb]:
                counters[lb] = 0
                current_indicator_sequence.loc[lb] = not current_indicator_sequence.loc[lb]
    return boolean_matrix    

# create oncoprint function

def oncoprint(
    color_matrix,
    marker_matrices = {},
    color_palette = {},
    sample_subset=None,
    feature_subset=None,
    row_colors=None,
    row_annot_width=0.05,
    row_annot_color_padding=0.01,
    row_annot_text_padding=(-0.5, 0.25),
    col_colors=None,
    col_annot_height=0.5,
    col_annot_color_padding=0.25,
    col_annot_text_padding=(-0.5, 0.25),
    show_column_annotation_label=True,
    figsize=(14, 8),
    show_frequency=[],
    label_frequency=True,
    row_sort='frequency',
    col_sort=['OncotreeLineage'],
    custom_column_order=None,
    alteration_sort=None,
    col_metadata=None,
    transparency=0.8,
    linewidth=1,
    circlesize=0.2,
    markerlinewidth=0.5,
):
    # subset the color matrix
    if sample_subset is None:
        sample_subset = color_matrix.index.tolist()
    if feature_subset is None:
        feature_subset = color_matrix.columns.tolist()
        
    # encode the values as integers which can be mapped to colors
    aligned_color_matrix = color_matrix.reindex(index=sample_subset, columns=feature_subset).T
    value_to_encoding = {v:i for i, v in enumerate(color_palette)}
    encoding_to_color = {i:color_palette[v] for v, i in value_to_encoding.items()}
    
    all_values = np.unique(aligned_color_matrix.values.ravel())
    aligned_color_matrix = aligned_color_matrix.replace(value_to_encoding)
    all_encodings = np.unique(aligned_color_matrix.values.ravel())
    sorted_colors = [encoding_to_color[i] for i in range(len(color_palette)) if i in all_encodings]              
                
    # assess the frequency of each feature
    frequency_fields = []
    combined_frequency = aligned_color_matrix.copy()
    combined_frequency.loc[:, :] = False
    for val in show_frequency:
        if val in all_values.tolist() or val in list(marker_matrices.keys()):
            frequency_fields.append(val)
    if len(frequency_fields) > 0:        
        for val in frequency_fields:
            if val in list(marker_matrices.keys()):
                combined_frequency = combined_frequency | (marker_matrices[val].reindex_like(aligned_color_matrix.T).T)
            elif val in all_values.tolist():
                combined_frequency = combined_frequency | (aligned_color_matrix == value_to_encoding[val])

        row_texts = (combined_frequency.mean(axis=1) * 100).apply(lambda x: f'{x:.1f}%')
    else:
        row_texts = None
        
    # assess the type overlap of each feature
    alteration_overlap_matrix = aligned_color_matrix.copy()
    alteration_overlap_matrix.loc[:, :] = 0
    all_alterations = all_values.tolist() + list(marker_matrices.keys())
    combination_component_encoding = pd.Series(
        [np.power(2, i) for i, x in enumerate(all_alterations)], index=all_alterations
    )
    indicator_map = generate_indicator_sequence_map(combination_component_encoding).reset_index().rename(
        {'index': 'indicator_sum'}, axis=1)
    if alteration_sort is not None:
        for i, alts in enumerate(alteration_sort):
            
            indicator_query = ' & '.join([f'`{x}`' for x in alts] + [f'~`{x}`' for x in all_alterations if x not in alts])
            indicator_match_idx = indicator_map.query(indicator_query).index.tolist()[0]
            indicator_map.loc[indicator_match_idx, 'order'] = i
        indicator_map['order'] = indicator_map['order'].fillna(len(alteration_sort))
    else:
        indicator_map['order'] = 0
    for val in all_alterations:
        if val in list(marker_matrices.keys()):
            alteration_overlap_matrix += (
                (marker_matrices[val].reindex_like(aligned_color_matrix.T).T).fillna(0).astype(int) * 
                combination_component_encoding[val]
            )
        elif val in all_values.tolist():
            alteration_overlap_matrix += (
                (aligned_color_matrix == value_to_encoding[val]).fillna(0).astype(int) * 
                combination_component_encoding[val]
            )
        
    # optionally sort rows by frequency
    if row_sort == 'frequency' and len(frequency_fields) > 0:
        row_order = combined_frequency.mean(axis=1).reset_index().sort_values([0, 'index'], ascending=[False, True])['index'].tolist()
        row_texts = row_texts.reindex(index=row_order)
        aligned_color_matrix = aligned_color_matrix.loc[row_order, :]
        combined_frequency = combined_frequency.reindex_like(aligned_color_matrix)
    # optionally sort columns by lineage, then presence of alterations
    if isinstance(col_sort, list) and len(col_sort) > 0:
        col_sort_order = col_sort + aligned_color_matrix.index.tolist()
        
        col_categoricals = dict()
        for c in col_sort:
            if c in col_metadata.columns.tolist():
                col_categoricals[c] = pd.Series(
                    pd.Categorical(
                        col_metadata[c].reindex(aligned_color_matrix.columns.tolist()),
                        categories=col_metadata[c].reindex(aligned_color_matrix.columns.tolist()).value_counts().index if custom_column_order is None else custom_column_order, 
                        ordered=True
                    ), index=aligned_color_matrix.columns.tolist(), name=c
                )
        
        if alteration_sort is not None:
            col_sort_df = pd.concat(
                list(col_categoricals.values()) + [alteration_overlap_matrix.T.replace(
                    indicator_map.set_index('indicator_sum')['order'].to_dict()
                )],axis=1
            )
            col_sort_ascending = [True for _ in col_categoricals] + [True for _ in alteration_overlap_matrix.index.tolist()]
        else:
            col_sort_df = pd.concat(
                list(col_categoricals.values()) + [combined_frequency.T],
                axis=1
            )
            col_sort_ascending = [True for _ in col_categoricals] + [False for _ in alteration_overlap_matrix.index.tolist()]

        col_order = col_sort_df.reset_index().sort_values(col_sort_order + ['index'], ascending=col_sort_ascending + [True])['index'].tolist()
        aligned_color_matrix = aligned_color_matrix.loc[:, col_order]
        combined_frequency = combined_frequency.reindex_like(aligned_color_matrix)
        alteration_overlap_matrix = alteration_overlap_matrix.reindex_like(aligned_color_matrix)
        
    # map the entries to plot coordinates
    feature_df = pd.DataFrame({
        'top': pd.Series(np.arange(len(feature_subset)), aligned_color_matrix.index.tolist()),
        'feature_center': pd.Series(np.arange(len(feature_subset)) + 0.5, aligned_color_matrix.index.tolist()),
        'bottom': pd.Series(np.arange(len(feature_subset)) + 1, aligned_color_matrix.index.tolist())
    })
    
    sample_df = pd.DataFrame({
        'left': pd.Series(np.arange(len(sample_subset)), aligned_color_matrix.columns.tolist()),
        'sample_center': pd.Series(np.arange(len(sample_subset)) + 0.5, aligned_color_matrix.columns.tolist()),
        'right': pd.Series(np.arange(len(sample_subset)) + 1, aligned_color_matrix.columns.tolist())
    })  
    
    # plot the base heatmap
    plt.figure(figsize=figsize)
    ax = sns.heatmap(aligned_color_matrix, cmap=sorted_colors, linewidths=linewidth, linecolor='black', cbar=False, xticklabels=True, yticklabels=True, alpha=transparency)
    
    # add column annotations
    if col_colors is not None:
        xlabel_order = [t.get_text() for t in ax.get_xticklabels()]
        if isinstance(col_colors, pd.Series):
            col_colors = pd.DataFrame(col_colors)
        col_colors = col_colors.reindex(index=xlabel_order).fillna('white')
        ax.tick_params(axis='x', which='major', top=False, bottom=False)
        for j, c in enumerate(col_colors.columns):
            for i, color in enumerate(col_colors[c]):
                ax.add_patch(plt.Rectangle(xy=(i, ((j+1) * -col_annot_height) - col_annot_color_padding), 
                                           width=1, height=col_annot_height, 
                                           facecolor=color, lw=linewidth, edgecolor='black',
                                           transform=ax.transData, clip_on=False))
            if show_column_annotation_label:
                ax.text(col_annot_text_padding[0], 0 - col_annot_height - col_annot_text_padding[1], c, transform=ax.transData)
            
    # move column labels to top
    ax.xaxis.set_label_position('top')
    ax.xaxis.tick_top()
    ax.set_xlabel('')
    ax.set_xticks([x-0.4 for x in ax.get_xticks()], aligned_color_matrix.columns.tolist(), rotation=60, ha='left')
    ax.tick_params(top=False)
    
    # move column labels if there are annotations present
    if col_colors is not None:
        for i, label in enumerate(ax.xaxis.get_majorticklabels()):
            label.set_transform(ax.transData)
            label.set_position((i, 0 - (2 * col_annot_height) - col_annot_text_padding[1]))
    
    # add marker-style annotations onto the base heatmap
    for marker_type, marker_matrix in marker_matrices.items():
        aligned_marker_matrix = marker_matrix.reindex_like(aligned_color_matrix.T).T
        long_marker_df = aligned_marker_matrix.melt(var_name='sample', ignore_index=False).reset_index().rename({'index': 'feature'}, axis=1)
        long_marker_df = long_marker_df.merge(feature_df, left_on='feature', right_index=True).merge(sample_df, left_on='sample', right_index=True)
        markers_to_plot = long_marker_df[long_marker_df['value'] == True]
        for _, r in markers_to_plot.iterrows():
            if marker_type == 'o':
                ax.add_patch(plt.Circle((r['sample_center'], r['feature_center']), circlesize, color='black'))
            elif marker_type == '\\':
                ax.plot([r['left'], r['right']], [r['top'], r['bottom']], color='black', linewidth=markerlinewidth) # backslash
            elif marker_type == '/':
                ax.plot([r['left'], r['right']], [r['bottom'], r['top']], color='black', linewidth=markerlinewidth) # forward slash
            elif marker_type == 'x':
                ax.plot([r['left'], r['right']], [r['top'], r['bottom']], color='black', linewidth=markerlinewidth) # backslash
                ax.plot([r['left'], r['right']], [r['bottom'], r['top']], color='black', linewidth=markerlinewidth) # forward slash
        
    # add frequency annotations
    if row_texts is not None and label_frequency:
        ax_right = ax.secondary_yaxis('right')
        row_positions = ax.get_yticks()
        ax_right.set_yticks(row_positions)
        ax_right.set_yticklabels(row_texts)
    
    return ax, pd.DataFrame().reindex_like(aligned_color_matrix)

# determine whether mutations are more frequent relative to a prior context

def categorize_mutations(
    mutation_matrix,
    group_subset,
    outgroup_subset,
    outgroup_max_models = 0,
    outgroup_max_fraction = 0.2,
    matched_annotations = None,
):
    mutation_matrix_subset = mutation_matrix.loc[group_subset, :]
    mutation_matrix_subset = mutation_matrix_subset.loc[:, mutation_matrix_subset.any().loc[lambda x: x].index.tolist()]
    outgroup_mutation_matrix_subset = mutation_matrix.reindex(index=outgroup_subset, columns=mutation_matrix_subset.columns.tolist())
    
    group_totals = mutation_matrix_subset.join(matched_annotations).groupby(matched_annotations.name).sum()
    outgroup_totals = outgroup_mutation_matrix_subset.join(matched_annotations).groupby(matched_annotations.name).sum()
    outgroup_n = outgroup_mutation_matrix_subset.join(matched_annotations).groupby(matched_annotations.name).count()
    outgroup_fractions = outgroup_totals.astype('Float64') / outgroup_n.astype('Float64')
    outgroup_matched_totals = outgroup_totals.reindex_like(group_totals)
    outgroup_matched_fractions = outgroup_fractions.reindex_like(group_totals)
    
    novel_by_matched_field = (
        (group_totals > 0) & 
        (outgroup_matched_totals <= outgroup_max_models) & 
        (outgroup_matched_fractions <= outgroup_max_fraction) #&
        # (outgroup_n > 0)
    )
    
    return novel_by_matched_field, group_totals, outgroup_matched_totals

def expand_annotation_to_models(collapsed_matrix, model_subset, model_metadata, expansion_column='OncotreeLineage'):
    group_to_models = model_metadata.loc[model_subset, expansion_column].groupby(
        model_metadata.loc[model_subset, expansion_column]
    ).groups
    output_matrix = []
    for c in collapsed_matrix.columns:
        for model in group_to_models[c]:
            output_matrix.append(collapsed_matrix.loc[:, c].rename(model))
    return pd.concat(output_matrix, axis=1)

def summarize_squares(alteration_matrix, shell, ingroup, outgroup, outgroup_max_n, model_metadata):
    alteration_locations = np.where(alteration_matrix.T.reindex_like(shell) == True)
    decorator_locations = np.where(
        expand_annotation_to_models(
            categorize_mutations(mutation_matrix=alteration_matrix, group_subset=ingroup, outgroup_subset=outgroup, outgroup_max_models = outgroup_max_n, matched_annotations=model_metadata.loc[:, 'OncotreeLineage'])[0].T, 
            model_subset=ingroup, model_metadata=model_metadata).reindex_like(shell).fillna(False) & 
        (alteration_matrix).T.reindex_like(shell).fillna(False)
    )
    
    alteration_table = pd.DataFrame({
        'Gene': shell.index[alteration_locations[0]],
        'GeneIndex': alteration_locations[0],
        'Model': shell.columns[alteration_locations[1]],
        'ModelIndex': alteration_locations[1],
    })
    decorator_table = pd.DataFrame({
        'Gene': shell.index[decorator_locations[0]],
        'GeneIndex': decorator_locations[0],
        'Model': shell.columns[decorator_locations[1]],
        'ModelIndex': decorator_locations[1],
        'Decorated': True
    })
    
    return alteration_table.merge(decorator_table, how='left').fillna(False)

def summarize_alteration_in_oncoprint(alteration_matrix, shell, ingroup, outgroup, rare_outgroup_max_n, model_metadata, alteration_type):
    rare_table = summarize_squares(alteration_matrix, shell, ingroup, outgroup, rare_outgroup_max_n, model_metadata).rename({'Decorated': 'RareInLineage'}, axis=1)
    novel_table = summarize_squares(alteration_matrix, shell, ingroup, outgroup, 0, model_metadata).rename({'Decorated': 'NovelInLineage'}, axis=1)
    return rare_table.merge(novel_table).assign(AlterationType=alteration_type).merge(model_metadata.loc[:, ['HasCRISPRData', 'OncotreeLineage']], left_on='Model', right_index=True)

def partition_models(omics_models_meta):
    next_gen_screened_models = omics_models_meta[
        omics_models_meta['IsNextGen'] &
        omics_models_meta['CRISPRScreenType'].isin(['3DO', '3DN', '2DN', '2DO']) &
        (~omics_models_meta['OncotreeLineage'].isin(['Other']))
    ].index.tolist()
    next_gen_screened_models_other = omics_models_meta[
        omics_models_meta['IsNextGen'] &
        omics_models_meta['CRISPRScreenType'].isin(['3DO', '3DN', '2DN', '2DO'])
    ].index.tolist()
    cell_line_screened_models = omics_models_meta[
        omics_models_meta['IsTraditional2D'] & 
        omics_models_meta['HasCRISPRData']
    ].index.tolist()
    
    next_gen_unscreened_models = omics_models_meta[
        omics_models_meta['IsNextGen'] & 
        (~omics_models_meta['HasCRISPRData'] | (omics_models_meta['CRISPRScreenType'] == '2DO')) & 
        omics_models_meta['HasCNData']].index.tolist()
    all_cell_line_models = omics_models_meta[
        omics_models_meta['IsTraditional2D'] & 
        omics_models_meta['HasCNData']
    ].index.tolist()
    
    return next_gen_screened_models, next_gen_unscreened_models, cell_line_screened_models, all_cell_line_models, next_gen_screened_models_other

def filter_genes(hotspot, next_gen_screened_models):
    geneset_table = load_file('geneset_table_updated.csv')
    og_genes = geneset_table[(geneset_table['Geneset'] == "Oncogenes")]['cds_gene_id'].tolist()
    tsg_genes = geneset_table[(geneset_table['Geneset'] == "Tumor suppressors")]['cds_gene_id'].tolist()
    
    # select only genes list
    genes_to_oncoprint = (
        set((hotspot.loc[next_gen_screened_models, :] > 0).sum().loc[lambda x: x > 0].index.tolist()) &
        (set(og_genes) | set(tsg_genes))
    )
    genes_to_oncoprint |= set(
        [search_gene(g, ) for g in ['CDH1', 'CDH2', 'RHOA', 'BRAF', 'ESR1']]
    )
    
    return genes_to_oncoprint


def create_large_oncoprint(cn_transform, hotspot, damaging, mutation_matrix, omics_models_meta, ingroup_models, outgroup_models, genes_to_oncoprint, lineage_order, figsize):

    # make the oncoprint
    ax, organoid_mtx_shell = oncoprint(
        color_matrix=cn_transform,
        col_metadata=omics_models_meta,
        color_palette=oncoprint_palette,
        marker_matrices={
            'o': hotspot > 0,
            '/': damaging > 0,
        },
        sample_subset=ingroup_models,
        feature_subset=list(genes_to_oncoprint),
        col_colors=omics_models_meta.loc[ingroup_models, 'OncotreeLineage'].replace(lineage_cmap),
        col_annot_height=0.5,
        col_annot_color_padding=0.25,
        col_annot_text_padding=(-6.5, -0.2),
        show_frequency=['Amplification', 'Deletion', '/', 'o'],
        alteration_sort=[['o', 'Amplification'], 
                         ['o', '/', 'Deletion'], 
                         ['o', '/', 'Neutral Copy'],
                         ['o', 'Neutral Copy'],
                         ['/', 'Deletion'],
                         ['/', 'Neutral Copy'],
                         ['Amplification'],
                         ['Deletion']],
        custom_column_order=lineage_order,
        figsize=figsize,
        transparency=0.7,
        markerlinewidth=1,
    )

    ax.set_yticks(ax.get_yticks(), [t.get_text().split(' ')[0] for t in ax.get_yticklabels()], fontsize=TICK_SIZE+1)
    ax.set_xticks(ax.get_xticks(), [t.get_text().split('-')[-1][-4:] for t in ax.get_xticklabels()], fontsize=TICK_SIZE+1)

    # identify the rare and novel alteration cases
    hotspot_table = summarize_alteration_in_oncoprint((hotspot > 0), organoid_mtx_shell, ingroup_models, outgroup_models, rare_outgroup_max_n=3, model_metadata=omics_models_meta, alteration_type='Hotspot')
    damaging_table = summarize_alteration_in_oncoprint((damaging > 0), organoid_mtx_shell, ingroup_models, outgroup_models, rare_outgroup_max_n=3, model_metadata=omics_models_meta, alteration_type='Damaging')
    amp_table = summarize_alteration_in_oncoprint((cn_transform == "Amplification"), organoid_mtx_shell, ingroup_models, outgroup_models, rare_outgroup_max_n=3, model_metadata=omics_models_meta, alteration_type='Amplification')
    del_table = summarize_alteration_in_oncoprint((cn_transform == "Deletion"), organoid_mtx_shell, ingroup_models, outgroup_models, rare_outgroup_max_n=3, model_metadata=omics_models_meta, alteration_type='Deletion')

    full_annotation_table = pd.concat([hotspot_table, damaging_table, amp_table, del_table], axis=0)
    full_annotation_table = full_annotation_table[~full_annotation_table['OncotreeLineage'].isin(['Other'])]
    rare_annotation_plotter_df = full_annotation_table[full_annotation_table['RareInLineage'] & ~full_annotation_table['NovelInLineage']].drop_duplicates(['GeneIndex', 'ModelIndex'])
    novel_annotation_plotter_df = full_annotation_table[full_annotation_table['NovelInLineage']].drop_duplicates(['GeneIndex', 'ModelIndex'])

    for i, j in zip(rare_annotation_plotter_df['GeneIndex'], rare_annotation_plotter_df['ModelIndex']):
        ax.add_patch(Rectangle((j, i), 1, 1, fill=False, edgecolor='green', lw=2.5, zorder=1e5-1, clip_on=False))

    for i, j in zip(novel_annotation_plotter_df['GeneIndex'], novel_annotation_plotter_df['ModelIndex']):
        ax.add_patch(Rectangle((j, i), 1, 1, fill=False, edgecolor='gold', lw=2.5, zorder=1e5, clip_on=False))

    # create the figure legend
    lineage_handles = ([Patch(facecolor='white', edgecolor='white')] + 
                       [Patch(facecolor=lineage_cmap[lin], edgecolor='black') for lin in lineage_order if lin not in ['Lung', 'Head and Neck']])
    lineage_labels = ['OncotreeLineage'] + [lin for lin in lineage_order if lin not in ['Lung', 'Head and Neck']]
    alteration_handles = [
        Patch(facecolor='white', edgecolor='white'),
        Patch(facecolor=plt.cm.coolwarm(256), edgecolor='black', linewidth=1, alpha=0.7),
        Patch(facecolor=plt.cm.coolwarm(0), edgecolor='black', linewidth=1, alpha=0.7),
        BackSlashRectangle(),
        (Patch(facecolor='white', edgecolor='black', linewidth=1), 
         Line2D([], [], color='black', marker='o', markersize=2, linewidth=0)),
        Patch(facecolor='white', edgecolor='gold', linewidth=1),
        Patch(facecolor='white', edgecolor='green', linewidth=1),
    ]
    alteration_labels = [
        'Genomic Alterations',
        'Amplification', 
        'Deletion', 
        'Damaging Mutation', 
        'Hotspot Mutation',
        'Novel Mutation\nin Lineage',
        'Rare Mutation\nin Lineage\n($\leq$ 3 Traditional)'
    ]
    ax.legend(
        handles=lineage_handles + alteration_handles, 
        labels=lineage_labels + alteration_labels,
        handlelength=1, handleheight=1,
        loc='lower left', bbox_to_anchor=(1.05, 0.5),
        handler_map=large_custom_handler_map
    )
    
    return full_annotation_table
    
# figure ED1a
    
def create_screened_oncoprint(cn_transform, hotspot, damaging, mutation_matrix, omics_models_meta, next_gen_screened_models, cell_line_screened_models, genes_to_oncoprint, lineage_order, filetype='pdf'):
    annotation_table = create_large_oncoprint(cn_transform, hotspot, damaging, mutation_matrix, omics_models_meta, ingroup_models=next_gen_screened_models, outgroup_models=cell_line_screened_models, genes_to_oncoprint=genes_to_oncoprint, lineage_order=lineage_order, figsize=(21, 6))

    if filetype == 'pdf':
        plt.savefig(os.path.join(FIGURE_DIR, 'Fig_ED1a_screened_next_gen_oncoplot.pdf'), bbox_inches='tight')
    if filetype == 'png':
        plt.savefig(os.path.join(FIGURE_DIR, 'Fig_ED1a_screened_next_gen_oncoplot.png'), bbox_inches='tight')
    annotation_table.to_csv(os.path.join(PROCESSED_DIR, 'screened_next_gen_oncoplot_annotation_summary.csv'), index=False)
    
# figure ED1a
    
def create_unscreened_oncoprint(cn_transform, hotspot, damaging, mutation_matrix, omics_models_meta, next_gen_unscreened_models, all_cell_line_models, genes_to_oncoprint, lineage_order, filetype='pdf'):
    annotation_table = create_large_oncoprint(cn_transform, hotspot, damaging, mutation_matrix, omics_models_meta, ingroup_models=next_gen_unscreened_models, outgroup_models=all_cell_line_models, genes_to_oncoprint=genes_to_oncoprint, lineage_order=lineage_order, figsize=(21, 6))

    if filetype == 'pdf':
        plt.savefig(os.path.join(FIGURE_DIR, 'Fig_ED1a_unscreened_next_gen_oncoplot.pdf'), bbox_inches='tight')
    elif filetype == 'png':
        plt.savefig(os.path.join(FIGURE_DIR, 'Fig_ED1a_unscreened_next_gen_oncoplot.png'), bbox_inches='tight')
    annotation_table.to_csv(os.path.join(PROCESSED_DIR, 'unscreened_next_gen_oncoplot_annotation_summary.csv'), index=False)

    
def get_alteration_shell(alteration_matrix, ingroup, outgroup, model_metadata, gene_subset, model_n=0):
    lineage_matrix = categorize_mutations(
        mutation_matrix=alteration_matrix, group_subset=ingroup, outgroup_subset=outgroup, outgroup_max_models=model_n, matched_annotations=model_metadata.loc[:, 'OncotreeLineage'],
    )[0].reindex(columns=gene_subset)
    lineage_matrix = lineage_matrix.loc[:, lineage_matrix.any()]
    
    lineage_matrix_expanded_to_model = expand_annotation_to_models(lineage_matrix.T, model_subset=ingroup, model_metadata=model_metadata)
    model_has_alteration = alteration_matrix.T.reindex_like(lineage_matrix_expanded_to_model)
    model_matrix = (model_has_alteration & lineage_matrix_expanded_to_model)
    model_matrix = model_matrix.loc[:, model_matrix.any()]
    return model_matrix, lineage_matrix
    
# figure 1c
    
def create_condensed_oncoprint(cn_transform, hotspot, damaging, mutation_matrix, omics_models_meta, genes_to_oncoprint, lineage_order):
    
    next_gen_screened_models, next_gen_unscreened_models, cell_line_screened_models, all_cell_line_models, next_gen_screened_models_other = partition_models(omics_models_meta)
    
    all_novel_mutations_by_model, all_novel_mutations_by_lineage = get_alteration_shell(mutation_matrix, next_gen_screened_models + next_gen_unscreened_models, all_cell_line_models, omics_models_meta, genes_to_oncoprint, model_n=0)
    all_novel_amplifications_by_model, all_novel_amplifications_by_lineage = get_alteration_shell((cn_transform == "Amplification"), next_gen_screened_models + next_gen_unscreened_models, all_cell_line_models, omics_models_meta, genes_to_oncoprint, model_n=0)
    all_novel_deletions_by_model, all_novel_deletions_by_lineage = get_alteration_shell((cn_transform == "Deletion"), next_gen_screened_models + next_gen_unscreened_models, all_cell_line_models, omics_models_meta, genes_to_oncoprint, model_n=0)
    
    condensed_novel_gene_list = list(
        set(all_novel_mutations_by_model.index.tolist() + all_novel_deletions_by_model.index.tolist() + all_novel_amplifications_by_model.index.tolist())
    )
    condensed_novel_model_list = list(
        set(all_novel_mutations_by_model.columns.tolist() + all_novel_deletions_by_model.columns.tolist() + all_novel_amplifications_by_model.columns.tolist())
    )

    # make the oncoprint
    ax, novel_mtx_shell = oncoprint(
        color_matrix=cn_transform,
        col_metadata=omics_models_meta,
        color_palette=oncoprint_palette,
        marker_matrices={
            'o': hotspot > 0,
            '/': damaging > 0,
        },
        sample_subset=condensed_novel_model_list,
        feature_subset=condensed_novel_gene_list,
        col_colors=omics_models_meta.loc[condensed_novel_model_list, ['CRISPRScreenType', 'OncotreeLineage']].replace(lineage_cmap).replace(
            {None: 'white', '2DO': 'white', '3DO': 'black', '2DN': 'black', '3DN': 'black'}),
        col_annot_height=0.8,
        col_annot_color_padding=0.5,
        col_annot_text_padding=(-20, -0.5),
        show_column_annotation_label=False,
        show_frequency=['Amplification', 'Deletion', '/', 'o'],
        label_frequency=False,
        alteration_sort=[['o', 'Amplification'], 
                         ['o', '/', 'Deletion'], 
                         ['o', '/', 'Neutral Copy'],
                         ['o', 'Neutral Copy'],
                         ['/', 'Deletion'],
                         ['/', 'Neutral Copy'],
                         ['Amplification'],
                         ['Deletion']],
        custom_column_order=lineage_order,
        figsize=(76 * mm, 48 * mm),
        transparency=0.7,
        linewidth=0.25,
        markerlinewidth=0.5,
        circlesize=0.1
    )

    # remove frequency ticks and add labels to the row annotations
    ax.set_xticklabels([])
    ax.set_yticks(ax.get_yticks(), [t.get_text().split(' ')[0] for t in ax.get_yticklabels()], fontsize=TICK_SIZE-1)
    ax.set_ylabel('')
    ax.set_title('')
    ax.text(-0.26, -1.9, 'Lineage', ha='right', va='center', fontdict={'fontsize': TICK_SIZE-0.5})
    ax.text(-0.26, -0.7, 'Screened', ha='right', va='center', fontdict={'fontsize': TICK_SIZE-0.5})

    # identify and annotate novel alterations
    hotspot_table = summarize_squares((hotspot > 0), novel_mtx_shell, condensed_novel_model_list, all_cell_line_models, outgroup_max_n=0, model_metadata=omics_models_meta).assign(AlterationType='Hotspot')
    damaging_table = summarize_squares((damaging > 0), novel_mtx_shell, condensed_novel_model_list, all_cell_line_models, outgroup_max_n=0, model_metadata=omics_models_meta).assign(AlterationType='Damaging')
    amp_table = summarize_squares((cn_transform == "Amplification"), novel_mtx_shell, condensed_novel_model_list, all_cell_line_models, outgroup_max_n=0, model_metadata=omics_models_meta).assign(AlterationType='Amplification')
    del_table = summarize_squares((cn_transform == "Deletion"), novel_mtx_shell, condensed_novel_model_list, all_cell_line_models, outgroup_max_n=0, model_metadata=omics_models_meta).assign(AlterationType='Deletion')
    full_annotation_table = pd.concat([hotspot_table, damaging_table, amp_table, del_table], axis=0).rename({'Decorated': 'NovelInLineage'}, axis=1).merge(omics_models_meta.loc[:, ['HasCRISPRData', 'OncotreeLineage']], left_on='Model', right_index=True)
    
    novel_annotation_plotter_df = full_annotation_table[full_annotation_table['NovelInLineage']].drop_duplicates(['GeneIndex', 'ModelIndex'])
    for i, j in zip(novel_annotation_plotter_df['GeneIndex'], novel_annotation_plotter_df['ModelIndex']):
        ax.add_patch(Rectangle((j, i), 1, 1, fill=False, edgecolor='gold', lw=0.75, zorder=1e5, clip_on=False))

    # create the figure legend
    alteration_handles = [
        Patch(facecolor=plt.cm.coolwarm(256), edgecolor='black', linewidth=0.5, alpha=0.7),
        Patch(facecolor=plt.cm.coolwarm(0), edgecolor='black', linewidth=0.5, alpha=0.7),
        BackSlashRectangle(linewidth=0.25),
        (Patch(facecolor='white', edgecolor='black', linewidth=0.5), 
         Line2D([], [], color='black', marker='o', markersize=2, linewidth=0)),
        Patch(facecolor='white', edgecolor='gold', linewidth=0.75),
    ]

    alteration_labels = [
        'Amplification', 
        'Deletion', 
        'Damaging mutation', 
        'Hotspot mutation',
        'Novel in lineage',
    ]

    ax.legend(
        handles=alteration_handles, 
        labels=alteration_labels,
        ncols=3, handlelength=1, handleheight=1, handletextpad=0.3, columnspacing=1,
        loc='lower left', bbox_to_anchor=(0.05, -0.21),
        fontsize=5,
        handler_map=custom_handler_map
    )
    plt.subplots_adjust(left=0.15, bottom=0.16, top=0.92, right=0.98)
    plt.savefig(os.path.join(FIGURE_DIR, 'Fig_1d_condensed_oncoplot.pdf'))
    full_annotation_table.to_csv(os.path.join(PROCESSED_DIR, 'condensed_oncoplot_annotation_summary.csv'), index=False)


def main():
    omics_models_meta, screen_metadata, damaging, hotspot, mutation_matrix, cn_transform = load_data()
    next_gen_screened_models, next_gen_unscreened_models, cell_line_screened_models, all_cell_line_models, next_gen_screened_models_other = partition_models(omics_models_meta)
    genes_to_oncoprint = filter_genes(hotspot, next_gen_screened_models)
    
    lineage_order = omics_models_meta.loc[next_gen_screened_models, 'OncotreeLineage'].value_counts().drop('CNS/Brain').index.tolist() + ['CNS/Brain']
    
    # figure 1d
    create_condensed_oncoprint(cn_transform, hotspot, damaging, mutation_matrix, omics_models_meta, genes_to_oncoprint, lineage_order)
    
    # figure ED1a
    lineage_order_with_other = omics_models_meta.loc[next_gen_screened_models, 'OncotreeLineage'].value_counts().drop('CNS/Brain').index.tolist() + ['CNS/Brain', 'Other']
    create_screened_oncoprint(cn_transform, hotspot, damaging, mutation_matrix, omics_models_meta, next_gen_screened_models_other, cell_line_screened_models, genes_to_oncoprint, lineage_order_with_other, filetype='pdf')

    # figure ED1a
    lineage_order_unscreened = omics_models_meta.loc[next_gen_unscreened_models, 'OncotreeLineage'].value_counts().drop('CNS/Brain').index.tolist() + ['CNS/Brain']
    create_unscreened_oncoprint(cn_transform, hotspot, damaging, mutation_matrix, omics_models_meta, next_gen_unscreened_models, all_cell_line_models, genes_to_oncoprint, lineage_order_unscreened, filetype='pdf')
    

if __name__ == "__main__":
    main()
    