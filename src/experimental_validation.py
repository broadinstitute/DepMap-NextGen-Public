import os
import pandas as pd
import numpy as np
from figure_utils import *
from constants import *
from data_utils import *

import seaborn as sns
import matplotlib.pyplot as plt

mpl.style.use(os.path.join(ASSETS_DIR, 'stylesheet.mplstyle'))

# figure 5h

def plot_viability_pointplots(viability_genes):
    viability_df = load_file('validation_normalized_viability.csv')
    viability_df['ModelCondition'] = viability_df['Model'] + ' ' + viability_df['Condition']

    viability_enumeration = {
        'KP3 Plastic': 0,
        'KP3 Dome-R10': 1,
        'KP3 Dome-OPAC': 2,
        'HS766T Plastic': 3,
        'HS766T Dome-R10': 4,
        'HS766T Dome-OPAC': 5
    }
    replicate_adj = {0: -0.2, 1: 0, 2: 0.2}

    viability_formatted_df = dict()
    for g in viability_df['Gene'].unique().tolist():
        viability_formatted_df[g] = viability_df[viability_df['Gene'] == g].copy()
        viability_formatted_df[g]['x'] = (
            viability_formatted_df[g]['ModelCondition'].replace(viability_enumeration).astype(float) + 
            viability_formatted_df[g]['Replicate'].replace(replicate_adj).astype(float)
        )
    viability_formatted_df = pd.concat(viability_formatted_df, ignore_index=True)
    
    fig, axs = plt.subplots(2, 2, figsize=(70 * mm, 55 * mm), gridspec_kw={'wspace':0.1, 'hspace':0.9}, sharey=True)
    axs = axs.flatten()

    plt.subplots_adjust(left=0.2, bottom=0.15, right=0.95, top=0.9)

    viability_genes = ['ITGB1', 'ITGAV', 'MSMO1', 'SQLE']

    for i, g in enumerate(viability_genes):
        viability_df_subset = viability_formatted_df[viability_formatted_df['Gene'] == g]
        sns.scatterplot(
            viability_df_subset,
            x='x', y='NormalizedViability',
            hue='Condition', palette=viability_palette,
            s=4,
            legend=False,
            ax = axs[i],
        )

        sns.pointplot(
            viability_df_subset,
            x='ModelCondition', y='NormalizedViability',
            color='black', 
            marker='_', markersize=4, markeredgewidth=0.5,
            linestyle='none', 
            errorbar='sd', err_kws={'linewidth': 0.25, 'color': 'black'}, capsize=0.25,
            legend=False,
            ax = axs[i],
        )

        axs[i].axhline(0, color='black', linestyle='dotted', linewidth=0.5)
        axs[i].axvline(2.5, color='gray', linestyle='solid', linewidth=0.75)
        axs[i].set_xticks(axs[i].get_xticks(), 
                          labels=[x.get_text().split(' ')[1] for x in axs[i].get_xticklabels()], 
                     rotation=30, horizontalalignment='right', fontsize=TICK_SIZE)
        axs[i].set_xlabel('')
        offset_x_labels(fig, axs[i], 20, 0)

        axs[i].set_title(g, fontsize=LABEL_SIZE, pad=8)
        axs[i].text(0.25, 1.05, 'KP3', fontdict={'fontsize': ANNOT_SIZE, 'color': 'black'}, horizontalalignment='center', transform=axs[i].transAxes)
        axs[i].text(0.75, 1.05, 'HS766T', fontdict={'fontsize': ANNOT_SIZE, 'color': 'black'}, horizontalalignment='center', transform=axs[i].transAxes)

    axs[2].set_ylabel('')
    axs[0].set_ylabel('Normalized viability', y=-0.5, fontdict={'fontsize': LABEL_SIZE})
    
    plt.savefig(os.path.join(FIGURE_DIR, 'viability_validation_pointplots.pdf'))

    
def main():
    # figure 5h
    plot_viability_pointplots((['ITGB1', 'ITGAV', 'MSMO1', 'SQLE']))

    
if __name__ == "__main__":
    main()   
    
