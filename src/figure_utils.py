import pandas as pd
import numpy as np

import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from matplotlib.patches import Patch, Rectangle
from matplotlib.lines import Line2D
from matplotlib.colors import ListedColormap, LinearSegmentedColormap, TwoSlopeNorm, Normalize, to_hex
from matplotlib.transforms import ScaledTranslation

from constants import *

TEXT_WHITE_BBOX_ARGS = {'facecolor': 'white', 'edgecolor': 'white', 'boxstyle':'round,pad=0.1'}

# Helper functions

def flatten_list(lst):
    return [x for sublist in lst for x in sublist]

def convert_categoricals_by_cmap(series, palette=dict(), null_value=0, null_color='tab:grey', palette_order=None):
    '''
    Replaces the entries in a series with an enumeration (integers), and additionally returns the mapping of values to integers and the ordered colormap
    '''
    all_values_in_series = series.unique().tolist()
    if palette_order is not None:
        ordered_values = [x for x in palette_order if x in all_values_in_series]
        ordered_values = ordered_values + sorted(list(set(all_values_in_series) - set(ordered_values)))
    else:
        ordered_values = all_values_in_series
    reverse_ordered_values = ordered_values[::-1]
    
    category_to_enum = dict()
    enum_to_color = dict()
    
    enum = null_value - 1
    for k in reverse_ordered_values:
        if k in palette:
            category_to_enum[k] = enum
            enum_to_color[enum] = palette[k]
            enum -= 1
        else:
            category_to_enum[k] = null_value
            enum_to_color[null_value] = null_color
    return series.replace(category_to_enum), category_to_enum, ListedColormap(
        [enum_to_color[k] for k in sorted(enum_to_color.keys())])

def prepare_gridspec(ax_map, height_ratios=None, width_ratios=None, figsize = (8, 10), 
                     inner_gs_dict=dict(), **outer_gs_kws):
    '''
    Creates a nested gridspec from a boolean matrix `ax_map` defining the subplots. This function only goes two layers deep in gridspecs.
    Named arguments to this function are applied to the outer gridspec, while inner gridspec arguments are supplied as a dictionary
    mapping coordinate positions in `ax_map` to dictionaries of GridSpecFromSubplotSpec keywords.
    '''
    fig = plt.figure(figsize=figsize)
    axis_matrix = ax_map.loc[ax_map.any(axis=1), ax_map.any()]
    nrows, ncols = axis_matrix.shape
    
    if height_ratios is None:
        height_ratios = [1 for _ in range(nrows)]
    elif width_ratios is None:
        width_ratios = [1 for _ in range(ncols)]
    
    outer = gridspec.GridSpec(nrows, ncols, width_ratios=width_ratios, height_ratios=height_ratios, **outer_gs_kws)
    ax_map = dict()
    for i in range(nrows):
        ax_map[i] = dict()
        for j in range(ncols):
            if axis_matrix.iloc[i, j]:
                if (i, j) in inner_gs_dict:
                    gs = gridspec.GridSpecFromSubplotSpec(subplot_spec=outer[i, j], **inner_gs_dict[(i, j)])
                    ax_map[i][j] = dict()
                    for mi in range(inner_gs_dict[(i, j)]['nrows']):
                        for mj in range(inner_gs_dict[(i, j)]['ncols']):
                            ax_map[i][j][(mi, mj)] = fig.add_subplot(gs[mi, mj])
                else:
                    gs = gridspec.GridSpecFromSubplotSpec(1, 1, subplot_spec=outer[i, j])
                    ax_map[i][j] = fig.add_subplot(gs[0])
    return fig, ax_map

def offset_x_labels(fig, ax, dx, dy):
    '''
    Adjusts x-axis labels horizontally in dpi units. Primarily used to fix alignment with rotated labels.
    '''
    offset = ScaledTranslation(dx / fig.dpi, dy / fig.dpi, fig.dpi_scale_trans)
    for label in ax.xaxis.get_majorticklabels():
        label.set_transform(label.get_transform() + offset)
        
def angle_label_rows(ax, text_ys, texts, anchor_ys, text_x = -0.5, **kwargs):
    '''
    Assigns a list of y-axis labels and creates angled tick marks toward them
    '''
    for text_y, label, anchor_y in zip(text_ys, texts, anchor_ys):
        ax.annotate(label, xy=(0, anchor_y), xycoords='data', 
                    xytext=(text_x, text_y), arrowprops={'arrowstyle': '-', 'color': 'black'}, annotation_clip=False, **kwargs)
        
def format_significance(value, intervals = {
    '***': pd.Interval(left=0, right=0.001, closed='both'),
    '**': pd.Interval(left=0.001, right=0.01, closed='right'),
    '*': pd.Interval(left=0.01, right=0.05, closed='right'),
    'n.s.': pd.Interval(left=0.05, right=1, closed='right')
}):
    '''
    Converts significance values to asterisk strings
    '''
    for label, interval in intervals.items():
        if value in interval:
            return label
    raise ValueError('input value does not fall within the provided intervals')
    
def approximate_fdr_in_p_units(df, fdr_threshold=0.05, pcolumn='pvalue', fdr_column='FDR'):
    least_sig_under = df[df[fdr_column] < fdr_threshold][pcolumn].max()
    most_sig_over = df[df[fdr_column] >= fdr_threshold][pcolumn].min()
    return (least_sig_under + most_sig_over) / 2

def search_annotations_index(texts, query_list):
    indices = []
    for i, q in enumerate(query_list):
        for j, t in enumerate(texts):
            if t.get_text() == q:
                indices.append(j)
        assert len(indices) > i, f"Text {q} was not found in annotation list"
    return indices

def manually_annotate(ax, text, df, xcol, ycol, xyoffset, xycoords='data', 
                      arrowprops=dict(arrowstyle='-', color='black', alpha=1, linewidth=0.25, shrinkA=1, shrinkB=1),
                      ha='center', va='center', color='black', weight='normal', 
                      bbox_kws=None):
    ax.annotate(
        text, (df.loc[text, xcol], df.loc[text, ycol]),
        xytext=(df.loc[text, xcol] + xyoffset[0], df.loc[text, ycol] + xyoffset[1]), xycoords=xycoords, 
        arrowprops=arrowprops, horizontalalignment=ha, verticalalignment=va,
        c=color, fontsize=5, fontweight=weight, bbox=bbox_kws
    )
    
def horizontal_annotate_significance_bar(ax, string, x, y_pair, text_offset=0.5, cap_length=0.5, linewidth=0.5):
    ax.plot([x, x], y_pair, color='black', linewidth=linewidth)
    ax.plot([x + cap_length, x], [y_pair[0], y_pair[0]], color='black', linewidth=linewidth)
    ax.plot([x + cap_length, x], [y_pair[1], y_pair[1]], color='black', linewidth=linewidth)
    if '*' in string: # adjust for vertical height of asterisks
        text_offset = text_offset / 3
    ax.text(
        x + text_offset, (y_pair[0] + y_pair[1]) / 2, string,
        horizontalalignment='center', verticalalignment='center',
        rotation=270,
        fontsize=ANNOT_SIZE, transform=ax.transData
    )
    
def dissect_legend(ax, expected_labels):
    '''
    Helper function to extract legend components by feature type.
    Arguments:
        ax: the matplotlib.pyplot.Axes object that holds the legend
        expected_labels: a list of strings that should match the header of each legend component. For example, if X is passed to the 'hue' argument of seaborn's
            scatterplot function, then it will automatically create a legend where X is the header for a set of colored patches following it
    Returns: a dictionary of header labels mapping to the tuple of [handles] and [labels] corresponding to it.
    '''
    if len(expected_labels) == 0:
        return dict()
    handles, labels = ax.get_legend_handles_labels()
    expected_label_breakpoints = dict()
    expected_label_breakpoints_r = dict()
    for i, lb in enumerate(labels):
        if lb in expected_labels:
            expected_label_breakpoints[lb] = i
            expected_label_breakpoints_r[i] = lb
    if len(expected_labels) == 1 and len(expected_label_breakpoints) == 0: 
        # if only one styling is used, matplotlib makes the legend header the title instead of a labeled blank patch, which cannot be recovered from (handles, labels)
        return {expected_labels[0]: (handles, labels)}
    breakpoints = list(expected_label_breakpoints.values())
    legend_components = dict()
    for j, bp in enumerate(breakpoints[:-1]):
        matching_label = expected_label_breakpoints_r[bp]
        legend_components[matching_label] = (handles[bp:breakpoints[j+1]], labels[bp:breakpoints[j+1]])
    legend_components[expected_label_breakpoints_r[breakpoints[-1]]] = (handles[breakpoints[-1]:], labels[breakpoints[-1]:])
    return legend_components

def bubble_plot(long_df, x_label, y_label, 
                xorder=None, yorder=None, 
                xcluster=None, ycluster=None,
                linkage_method='single',
                hue_label=None, hue_min=np.nan, hue_max=np.nan, hue_center=None, palette='RdBu_r',
                size_label=None, size_min=np.nan, size_max=np.nan, 
                style_label=None,
                border_label=None, linewidth=1,
                figsize=(8, 6), legend_loc=(1.15, 0.6), width_ratios=[29, 1], gridspec_kws={}, **plt_kws):
    '''
    Creates a variation of a heatmap, where each gridpoint can vary by hue, size, style, and border. 
    Arguments:
        long_df: a long-form pandas.DataFrame containing the data to be plotted
        x_label, y_label: the column name corresponding to the x-axis, y-axis
        xorder, yorder: lists of entries in the x_label, y_label columns. If provided, the axes will be arranged in the given order. Otherwise, the row and column
            orders will be as they appear in the input 'long_df'
        xcluster, ycluster: the column names corresponding to the values to be used for clustering. If provided, converts 'long_df' into a matrix with x_label on
            the rows, y_label on the columns, and xcluster (ycluster) as the values, then performs hierarchical clustering along the axis. Accepts "hue", "size",
            and "style", in which case the column name will be inferred from the corresponding arguments.
        linkage_method: the hierarchical clustering method to be used. See 'scipy.cluster.hierarchy.linkage' for options.
        hue_label: the column name corresponding to the color variable
        hue_min, hue_max: the bounds for the hue values in data units. Values in the hue column will be clipped to fit the bounds specified
        hue_center: if provided, creates a diverging colorscale from the specified center
        palette: the colormap to use
        size_label: the column name corresponding to the size variable
        size_min, size_max: the bounds for the size variable in data units. Not to be confused with the plotting argument "sizes", which controls the point size
        style_label: the column name corresponding to the style variable
        border_label: the column name corresponding to the border variable. MUST be a boolean column, where True entries will receive a border
        linewidth: the border width around points
        figsize: the size of the matplotlib.pyplot.Figure to create
        legend_loc: passed to 'bbox_to_anchor' in 'matplotlib.pyplot.legend' to manually orient the legend. In units of axes transform, with reference to the lower
            left of the figure.
    Returns:
        the matplotlib.pyplot.Axes object for the figure and the matplotlib.colorbar.Colorbar object, if hue specified
    '''

    from matplotlib.colors import TwoSlopeNorm, Normalize
    from scipy.cluster.hierarchy import linkage, dendrogram
    from warnings import warn

    long_df = long_df.copy()

    # preprocess the coloring
    if hue_label is not None:
        if hue_min is None or np.isnan(hue_min):
            hue_min = long_df[hue_label].min()
        if hue_max is None or np.isnan(hue_max):
            hue_max = long_df[hue_label].max()
        if hue_center is not None: # make a divergent palette if the center is specified
            color_norm = TwoSlopeNorm(
                vcenter=hue_center, 
                vmin=hue_min,
                vmax=hue_max,
            )
        else: # else make a continuous palette
            color_norm = Normalize(
                vmin=hue_min,
                vmax=hue_max
            )
    else:
        color_norm = None

    # preprocess the sizes
    if size_label is not None:
        if size_min is None or np.isnan(size_min):
            size_min = long_df[size_label].min()
        else:
            long_df[size_label] = long_df[size_label].clip(lower=size_min)
        if size_max is None or np.isnan(size_max):
            size_max = long_df[size_label].max()
        else:
            long_df[size_label] = long_df[size_label].clip(upper=size_max)
        size_norm = Normalize(
            vmin=size_min,
            vmax=size_max
        )
    else:
        size_norm=None

    # fix the x-axis and y-axis ordering by mapping the unique elements of their correponding columns to coordinates
    if xorder is None: # if unspecified, use the order found in the input dataframe or cluster the x-axis
        if xcluster is None:
            xorder = long_df[x_label].tolist()
        else:
            if xcluster == 'hue' and hue_label is not None:
                xcluster = hue_label
            elif xcluster == 'size' and size_label is not None:
                xcluster = size_label
            if ycluster == 'style' and style_label is not None:
                style_enumeration = pd.Series({v:k for k,v in enumerate(pd.unique(long_df[style_label].unique().tolist()))})
                xmtx_data = long_df.pivot_table(index=x_label, columns=y_label, values=style_label).replace(style_enumeration)
            else:
                xmtx_data = long_df.pivot_table(index=x_label, columns=y_label, values=xcluster)
            if xmtx_data.isna().sum().sum() > 0:
                warn('Clustering matrix for columns contains NaN values. Filling NaNs with 0 for clustering only...')
                xmtx_data = xmtx_data.fillna(0)
            xdata_linkage = linkage(xmtx_data, method=linkage_method)
            xcluster_order = dendrogram(xdata_linkage, no_plot=True)['leaves']
            xorder = xmtx_data.iloc[xcluster_order, :].index.tolist()
    xorder = pd.Series({v:k for k,v in enumerate(pd.unique(xorder).tolist())})

    if yorder is None: # if unspecified, use the order found in the input dataframe or cluster the y-axis
        if ycluster is None:
            yorder = long_df[y_label].tolist()
        else:
            if ycluster == 'hue' and hue_label is not None:
                ycluster = hue_label
            elif ycluster == 'size' and size_label is not None:
                ycluster = size_label
            if ycluster == 'style' and style_label is not None:
                style_enumeration = pd.Series({v:k for k,v in enumerate(sorted(long_df[style_label].unique().tolist()))})
                ymtx_data = long_df.pivot_table(index=y_label, columns=x_label, values=style_label).replace(style_enumeration)
            else:
                ymtx_data = long_df.pivot_table(index=y_label, columns=x_label, values=ycluster)
            if ymtx_data.isna().sum().sum() > 0:
                warn('Clustering matrix for rows contains NaN values. Filling NaNs with 0 for clustering only...')
                ymtx_data = ymtx_data.fillna(0)
            ydata_linkage = linkage(ymtx_data, method=linkage_method)
            ycluster_order = dendrogram(ydata_linkage, no_plot=True)['leaves']
            yorder = ymtx_data.iloc[ycluster_order, :].index.tolist()
    yorder = pd.Series({v:k for k,v in enumerate(pd.unique(yorder).tolist())})

    # filter dataframe to only elements that are explicitly ordered
    long_df = long_df[long_df[x_label].isin(xorder.index.tolist()) & long_df[y_label].isin(yorder.index.tolist())]
    if np.any(np.array([len(x) for x in long_df.groupby([x_label, y_label]).groups.values()]) > 1):
        warn(f'Dataframe contains duplicate entries for ({x_label}, {y_label}) pairs. They will be overplotted in order of the inputs...')

    # map the coordinates
    long_df['x_coordinate'] = long_df[x_label].replace(xorder)
    long_df['y_coordinate'] = long_df[y_label].replace(yorder)

    # make the plot
    fig, axs = plt.subplots(1, 2, figsize=figsize, width_ratios=width_ratios, gridspec_kw=gridspec_kws)
    
    if border_label is not None:
        sns.scatterplot(
            long_df[~long_df[border_label]],
            x='x_coordinate', y='y_coordinate',
            hue=hue_label, hue_norm=color_norm, palette=palette,
            size=size_label, size_norm=size_norm,
            style=style_label,
            linewidth=0,
            ax=axs[0],
            **plt_kws
        )
        sns.scatterplot(
            long_df[long_df[border_label]],
            x='x_coordinate', y='y_coordinate',
            hue=hue_label, hue_norm=color_norm, palette=palette,
            size=size_label, size_norm=size_norm,
            style=style_label,
            linewidth=linewidth,
            ax=axs[0],
            **plt_kws
        )
    else:
        sns.scatterplot(
            long_df,
            x='x_coordinate', y='y_coordinate',
            hue=hue_label, hue_norm=color_norm, palette=palette,
            size=size_label, size_norm=size_norm,
            style=style_label,
            linewidth=linewidth,
            ax=axs[0],
            **plt_kws
        )


    # add colorbar
    if hue_label is not None:
        sm = plt.cm.ScalarMappable(cmap=palette, norm=color_norm)
        sm.set_array(long_df[hue_label].dropna())

        cbar = axs[1].figure.colorbar(sm, cax=axs[1])
        cbar.ax.set_ylabel(hue_label, rotation=270, labelpad=15)
    else:
        cbar = None
        axs[1].set_visible(False)

    # formatting
    axs[0].grid(True)
    axs[0].set_axisbelow(True)
    axs[0].set_xticks(xorder.values.tolist(), xorder.index.tolist(), rotation=90);
    axs[0].set_yticks(yorder.values.tolist(), yorder.index.tolist());
    axs[0].set_xlim(-1, len(xorder))
    axs[0].set_ylim(-1, len(yorder))

    # remake the legend by removing the color component and adding back the size and style components, if present
    legend_elements = [z for z in [hue_label, size_label, style_label] if z is not None]
    legend_components = dissect_legend(axs[0], legend_elements)
    handles_to_show = []
    labels_to_show = []
    if size_label in legend_elements:
        handles_to_show.extend(legend_components[size_label][0])
        labels_to_show.extend(legend_components[size_label][1])
    if style_label in legend_elements:
        h_to_add = legend_components[style_label][0]
        for h in h_to_add:
            h.set_markersize(10)
        handles_to_show.extend(legend_components[style_label][0])
        labels_to_show.extend(legend_components[style_label][1])
    if len(legend_elements) == 1: # case where one styling is used, need to set the legend title manually because it is not a labeled blank patch
        axs[0].legend(title=legend_elements[0], 
                    handles=handles_to_show, labels=labels_to_show, 
                    loc='lower left', bbox_to_anchor=legend_loc)
    elif len(legend_components) > 0:
        axs[0].legend(handles=handles_to_show, labels=labels_to_show, 
                    loc='lower left', bbox_to_anchor=legend_loc)

    return axs, cbar

def make_gene_effect_distribution_plot(genes, mtx=None,
                                       x='Gene', y='Gene Effect',
                                       sample_subset=None,
                                       grouping_annotation=None, figsize=None, 
                                       point_hue=None, point_palette=None, point_alpha=None,
                                       box_hue=None, box_palette=None, box_alpha=None,
                                       dodge=True, hue_order=None, plot_type='boxstrip', 
                                       box_kwargs={}, point_kwargs={}, ax=None):
    if sample_subset is None:
        sample_subset = mtx.index.tolist()
    if type(genes) == str:
        genes = [genes]
    full_genes = genes
    to_plot = pd.DataFrame(mtx.loc[sample_subset, full_genes]).melt(
        ignore_index=False, var_name='Gene', value_name='Gene Effect'
    ).join(grouping_annotation).dropna()
    
    if ax is None:
        if figsize is None:
            figsize=(6, 4)
        plt.figure(figsize=figsize)
        ax = plt.gca()
    
    if x is None:
        x = 'Gene'
    if y is None:
        y = 'Gene Effect'
    if point_hue is None:
        point_hue = grouping_annotation.name
    if box_hue is None:
        box_hue = point_hue
    
    if plot_type == 'boxen':
        sns.boxenplot(to_plot, x=x, y=y, hue=box_hue, palette=box_palette, hue_order=hue_order, ax=ax, dodge=dodge, alpha=box_alpha, **box_kwargs)
    elif plot_type == 'strip':
        sns.stripplot(to_plot, x=x, y=y, hue=point_hue, palette=point_palette, hue_order=hue_order, ax=ax, dodge=dodge, alpha=point_alpha, **point_kwargs)
    elif plot_type == 'swarm':
        sns.swarmplot(to_plot, x=x, y=y, hue=point_hue, palette=point_palette, hue_order=hue_order, ax=ax, dodge=dodge, alpha=point_alpha, **point_kwargs)
    elif plot_type == 'box':
        sns.boxplot(to_plot, x=x, y=y, hue=box_hue, palette=box_palette, hue_order=hue_order, ax=ax, dodge=dodge, alpha=box_alpha, **box_kwargs)
    elif 'strip' in plot_type and 'box' in plot_type:
        sns.boxplot(to_plot, x=x, y=y, hue=box_hue, palette=box_palette, hue_order=hue_order, ax=ax, dodge=dodge, alpha=box_alpha,
                    showfliers=False, **box_kwargs)
        sns.stripplot(to_plot, x=x, y=y, hue=point_hue, palette=point_palette, hue_order=hue_order, ax=ax, dodge=dodge, alpha=point_alpha,
                      **point_kwargs)
    elif 'violin' in plot_type and 'box' in plot_type:
        sns.violinplot(to_plot, x=x, y=y, hue=box_hue, palette=box_palette, hue_order=hue_order, ax=ax, dodge=dodge, alpha=box_alpha,
                       **box_kwargs)
        sns.stripplot(to_plot, x=x, y=y, hue=point_hue, palette=point_palette, hue_order=hue_order, ax=ax, dodge=dodge, alpha=point_alpha,
                      **point_kwargs)
        
    return ax