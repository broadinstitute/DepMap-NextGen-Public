import os
import pandas as pd
import numpy as np
from tqdm import tqdm

import scipy.stats
import statsmodels.api as sm
import statsmodels.formula.api as smf

from data_utils import *
from constants import *
from gene_utils import *

from chronos import Chronos, normalize_readcounts, calculate_fold_change
from sklearn.metrics import roc_curve, auc

import gseapy

SEED = 1
CHRONOS_SEED = 41

### Perform gene classifications on dependency scores ###

def tda_classify_common_essentials(gene_effect_matrix, quantile=0.9):
    print('Identifying common essentials...')
    rank_data = gene_effect_matrix.rank(axis=1, method='min', na_option='bottom') # Rank the data across rows
    rank_data = rank_data.div(gene_effect_matrix.notna().sum(axis=1), axis=0) # Divide by the number of non-NA values for each row to normalize
    percentile = rank_data.quantile(q=quantile) # Compute the 90th percentile for each column in the DataFrame
    
    # Calculate the density function of the 90th percentile values
    kde = scipy.stats.gaussian_kde(percentile.dropna(), bw_method=0.3)
    test_points = np.arange(0.1, 0.9, 0.001)
    threshold = test_points[np.argmin(kde(test_points))]
    ce_list = pd.DataFrame(
        {"Row.name": percentile.index, "CE_percentile": percentile.values}
    )
    
    print(f'Common essentials are in the {threshold * 100:.1f}% most depleted genes')
    
    # Add a boolean column to the DataFrame that indicates whether the CE_percentile value is <= to the threshold
    ce_list["Common_Essential"] = ce_list["CE_percentile"] <= threshold
    ce_list = ce_list.set_index('Row.name')
    print(f"Weakest common essential has mean {gene_effect_matrix.loc[:, ce_list.query('Common_Essential').index].mean().max():.2f}")
    
    return ce_list

def tda_classify_strongly_selective(gene_effect_matrix, probability_dependency_matrix, probability_threshold=0.5, selectivity_threshold=-0.86):
    print('Identifying strongly selective genes...')
    mean = gene_effect_matrix.mean(axis=0)  # Calculate the mean of each column
    median = gene_effect_matrix.median(axis=0)  # Calculate the median of each column
    var = gene_effect_matrix.var(axis=0)  # Calculate the variance of each column
    skewness = gene_effect_matrix.skew(axis=0)  # Calculate the skewness of each column
    kurtosis = gene_effect_matrix.kurtosis(axis=0)  # Calculate the skewness of each column
    n = (1 - kurtosis.isna()).sum()
    bimodality = (skewness ** 2 + 1) / (kurtosis + 3 * (n - 1) ** 2 / ((n - 2) * (n - 3)))

    moments = pd.DataFrame(
        {
            "Mean": mean,
            "Median": median,
            "Variance": var,
            "Skewness": skewness,
            "Kurtosis": kurtosis,
            "Bimodality": bimodality,
        }
    )
    
    dep_lines = (probability_dependency_matrix > probability_threshold).apply(sum)
    lines_with_data = (~pd.isna(probability_dependency_matrix)).apply(sum)
    merged = pd.DataFrame(dict(dep_lines=dep_lines, lines_with_data=lines_with_data))
    merged = moments.join(merged)
    
    merged['skew_x_kurtosis'] = merged["Skewness"] * merged["Kurtosis"]
    merged["is_strongly_selective"] = (merged['skew_x_kurtosis'] < selectivity_threshold) & (merged["dep_lines"] > 0)
    
    return merged

def tda_classify_high_variance(gene_dependency_matrix, expr_matrix):
    print('Identifying high variance genes...')
    common_models = sorted(list(set(gene_dependency_matrix.index) & set(expr_matrix.index)))
    common_genes = sorted(list(set(gene_dependency_matrix.columns) & set(expr_matrix.columns)))
    
    dep_mtx = gene_dependency_matrix.loc[common_models, :]
    expr_mtx = expr_matrix.loc[common_models, common_genes]
    
    # identify the modes of the gene expression distribution
    kde = scipy.stats.gaussian_kde(expr_mtx.melt()['value'].dropna())
    test_points = np.arange(1, 4, .005)
    expr_cutoff = test_points[np.argmin(kde(test_points))]
    
    nonexp_frac = (expr_mtx < expr_cutoff).mean()
    nonexp_genes = nonexp_frac.loc[lambda x: x == 1].index.tolist()
    
    print(f'{len(nonexp_genes)} genes found as unexpressed (TPM < {expr_cutoff:.2f}) across all models')

    var_df = pd.DataFrame(dep_mtx.var().rename('prob_variance'))
    var_df['unexpressed'] = var_df.index.isin(nonexp_genes)
    
    # use the unexpressed genes to define a variance threshold
    var_thresh = var_df.query('unexpressed')['prob_variance'].quantile(0.99)
    var_df['high_variance'] = (var_df['prob_variance'] > var_thresh) & (~var_df['unexpressed'])
    
    return var_df

def run_tda_classification(gene_effect_matrix, gene_dependency_matrix, expr_matrix, model_subsets = dict(), 
                           probability_threshold = 0.5, selectivity_threshold=-0.86, essentiality_quantile=0.9):
    results_dict = dict()
    for subset_name, ms in tqdm(model_subsets.items()):
        print(f'Using {len(ms)} models...')
        print(f'Classifying {subset_name}...')
        tda_essential = tda_classify_common_essentials(gene_effect_matrix.loc[ms, :], quantile=essentiality_quantile)
        tda_selective = tda_classify_strongly_selective(gene_effect_matrix.loc[ms, :], gene_dependency_matrix.loc[ms, :], 
                                                        probability_threshold=probability_threshold, selectivity_threshold=selectivity_threshold)
        tda_high_variance = tda_classify_high_variance(gene_dependency_matrix.loc[ms, :], expr_matrix)
        results_dict[subset_name] = tda_essential.join(tda_selective).join(tda_high_variance)
    return results_dict

def reclassify_tda_selective(results_df, min_frac_dep=None, min_n_dep=1, threshold=-0.86):
    results_copy = results_df.copy()
    results_copy['frac_dep'] = results_copy['dep_lines'] / results_copy['lines_with_data']
    if min_frac_dep is not None:
        min_deps_mask = results_copy['frac_dep'] >= min_frac_dep
    elif min_n_dep is not None:
        min_deps_mask = results_copy['dep_lines'] >= min_n_dep
    reclassified = ((results_copy['skew_x_kurtosis'] < threshold) & min_deps_mask).loc[lambda x: x].index.tolist()
    results_copy['selective'] = results_copy.index.isin(reclassified)
    return results_copy

def partition_genes_by_tda_results(tda_results_df, selectivity_kws={}, paness_threshold=None,
                                   order = ['Non-dependency', 'Weakly Selective', 'Pan-dependency', 
                                            'High Variance', 'Strongly Selective']):
    new_calls_full = reclassify_tda_selective(tda_results_df, **selectivity_kws)
    unscreened = new_calls_full.query('lines_with_data == 0').index.tolist()
    
    new_calls = new_calls_full.drop(unscreened, axis=0)
    class_to_genes = dict()
    class_to_genes['Strongly Selective'] = new_calls.query('selective').index.tolist()
    class_to_genes['High Variance'] = new_calls.query('high_variance').index.tolist()
    if paness_threshold is not None:
        class_to_genes['Pan-dependency'] = new_calls.query(f'CE_percentile <= {paness_threshold}').index.tolist()
    else:
        class_to_genes['Pan-dependency'] = new_calls.query('Common_Essential').index.tolist()
    class_to_genes['Weakly Selective'] = new_calls.query('dep_lines > 0').index.tolist()
    class_to_genes['Non-dependency'] = new_calls.query('dep_lines == 0').index.tolist()
    
    tda_class = pd.Series('Unscreened', index=new_calls_full.index.tolist())
    for tc in order:
        tda_class.loc[class_to_genes[tc]] = tc
    
    return pd.DataFrame({
        'Class': tda_class,
        'Variance': new_calls_full['Variance'],
        'Selectivity': new_calls_full['skew_x_kurtosis'],
        'DepLines': new_calls_full['dep_lines'],
        '90%ileLineRank': new_calls_full['CE_percentile']
    })

def call_gene_classes():
    # load data
    screen_gene_effect = load_file('screen_gene_effect.csv', index_col=0)
    screen_gene_dependency = load_file('screen_gene_dependency.csv', index_col=0)
    screen_metadata = load_file('screen_metadata.csv', index_col=0)
    omics_models_meta = load_file('model_metadata.csv', index_col=0)
    expr = load_full_matrix('expression')
    screen_expr = expand_model_matrix_to_screens(expr, screen_metadata.reset_index().set_index('ModelID')['ScreenID'])
    
    # define model sets
    organoid_screens = screen_metadata[
        screen_metadata['ScreenType'].isin(['2DO', '3DO']) # define new set class
    ]
    organoid_screens = organoid_screens.loc[
        organoid_screens.index.isin(screen_gene_effect.index)
    ]
    next_gen_screens = screen_metadata[
        screen_metadata['ScreenType'].isin(['2DO', '3DO', '3DN', '2DN'])
    ]
    next_gen_screens = next_gen_screens.loc[
        next_gen_screens.index.isin(screen_gene_effect.index)
    ]
    next_gen_lineages = next_gen_screens.value_counts('OncotreeLineage').index.tolist()
    
    # run the gene classification
    screen_tda_classifications = run_tda_classification(
        screen_gene_effect, screen_gene_dependency, screen_expr,
        model_subsets={
            'Organoids': organoid_screens.index.tolist(),
            'NextGen': next_gen_screens.index.tolist()
        }
    )
    
    # compile results
    tda_summary = pd.concat([
        partition_genes_by_tda_results(
            dataset[grp], selectivity_kws={'min_n_dep': 3},
            order=['Non-dependency', 'Weakly Selective', 'High Variance', 'Strongly Selective', 'Pan-dependency']
        ).add_prefix(grp)
        for dataset, grp in zip(
            [screen_tda_classifications, screen_tda_classifications],
            ['NextGen', 'Organoids']
        )
    ], axis=1)
    
    tda_summary.to_csv(os.path.join(PROCESSED_DIR, 'tda_dependency_classification.csv'))
    
### Compute metaprogram scores and celligner metadata ###

# normalize expression data according to Tirosh, 2016
def z_score_by_expression_bins(expression_matrix, nbin=25, n_sample=100, seed=None):
    # compute gene means and bin them
    mean_expr_profile = expression_matrix.mean()
    genes_to_bins = pd.cut(mean_expr_profile, bins=nbin)
    
    if seed is not None:
        np.random.seed(seed)
        
    # identify random control sets of genes in each bin
    control_genes_by_bin = genes_to_bins.reset_index().rename({'index': 'gene', 0: 'bin'}, axis=1).groupby(
        'bin'
    ).apply(lambda x: x.sample(n=100) if len(x) >= 100 else x.sample(frac=1)).set_index('gene')['bin']
    genes_in_use = pd.concat([mean_expr_profile.rename('mean'), genes_to_bins.rename('bin')], axis=1)
    genes_in_use['control'] = genes_in_use.index.isin(control_genes_by_bin.index.tolist())
    
    control_means = expression_matrix.groupby(control_genes_by_bin, axis=1).mean()
    
    # normalize gene expression values against control means in respective expression bins
    aligned_control_means = control_means.loc[:, genes_to_bins.loc[expression_matrix.columns.tolist()].values]
    normalized_expression_matrix = pd.DataFrame(
        expression_matrix.values - aligned_control_means.values,
        index=expression_matrix.index.tolist(), columns=expression_matrix.columns.tolist()
    )
    
    return normalized_expression_matrix, genes_in_use

def compute_geneset_scores():
    # Celligner data - import and map the gene ids to the right format
    celligner_coordinates = load_file('celligner_coordinates_w_hcmi_filtered.csv', index_col=0)
    tumor_expression = load_file('celligner_corrected_expr_w_hcmi_filtered.csv.gz', compression='gzip', index_col=0)
    omics_models_meta = load_file('model_metadata.csv', index_col=0)
    expr = load_full_matrix('expression')
    geneset_table = load_file('geneset_table_updated.csv')
    
    # Gavish
    mmps_dict = geneset_table[geneset_table['Source'] == "Gavish"].groupby('Geneset').apply(lambda x: list(x['cds_gene_id'])).to_dict()
    
    # Moffitt
    basal_genes = geneset_table[(geneset_table['Source'] == 'Moffitt') & (geneset_table['Geneset'] == "Basal")]['cds_gene_id'].tolist()
    classical_genes = geneset_table[(geneset_table['Source'] == 'Moffitt') & (geneset_table['Geneset'] == "Classical")]['cds_gene_id'].tolist()
    
    all_celligner_samples = celligner_coordinates[
        (celligner_coordinates['type'].isin(['TCGA+ tumor', 'MET500 tumor', 'Novartis_PDX'])) |
        (celligner_coordinates.index.isin(
            omics_models_meta[omics_models_meta['HasExprData']].index.tolist()
        ))
    ]
    
    # normalize expression datasets separately (DepMap 2D+3D, TCGA)
    depmap_expr = expr.reindex(index=omics_models_meta[omics_models_meta['HasExprData']].index.tolist()).dropna(how='all')
    tcga_expr = tumor_expression.reindex(
        index=all_celligner_samples[all_celligner_samples['type'] == 'TCGA+ tumor'].index.tolist()
    )
    
    np.random.seed(SEED) #HERE
    depmap_zscore_expr, depmap_zscore_controls = z_score_by_expression_bins(depmap_expr, nbin=25, n_sample=100, seed=None)
    tcga_zscore_expr, tcga_zscore_controls = z_score_by_expression_bins(tcga_expr, nbin=25, n_sample=100, seed=None)

    # compute scores for metaprograms
    depmap_mmp_scores = pd.concat(
        [depmap_zscore_expr.reindex(columns=mmps_dict[mmp]).mean(axis=1).rename(mmp) for mmp in mmps_dict], axis=1
    )
    tcga_mmp_scores = pd.concat(
        [tcga_zscore_expr.reindex(columns=mmps_dict[mmp]).mean(axis=1).rename(mmp) for mmp in mmps_dict], axis=1
    )
    
    # compute basal and classical scores
    moffitt_scores = pd.concat([
        pd.concat([
            depmap_zscore_expr.reindex(columns=basal_genes).mean(axis=1),
            tcga_zscore_expr.reindex(columns=basal_genes).mean(axis=1),
        ], axis=0).rename('Moffitt Basal'),
        pd.concat([
            depmap_zscore_expr.reindex(columns=classical_genes).mean(axis=1),
            tcga_zscore_expr.reindex(columns=classical_genes).mean(axis=1),
        ], axis=0).rename('Moffitt Classical')
    ], axis=1)
    
    # merge Gavish scores into a matrix
    all_mmp_scores = pd.concat([depmap_mmp_scores, tcga_mmp_scores], axis=0)
    # merge Moffitt scores into the same matrix
    all_geneset_zscores = all_mmp_scores.join(moffitt_scores, how='outer')
    

    # TODO: REMOVE HACKY SOLN
    # keep rows that start with ACH- 
    # replace other rows with previous celligner_geneset_scores from 2025 analytical results
    prev_tcga_scores = pd.read_csv(
        "tmp-files/tcga_geneset_scores.csv", index_col= 0 
    )
    tmp = prev_tcga_scores.combine_first(all_geneset_zscores)

    tmp.to_csv(os.path.join(PROCESSED_DIR, 'celligner_geneset_scores.csv'))

    # all_geneset_zscores.to_csv(os.path.join(PROCESSED_DIR, 'celligner_geneset_scores.csv'))


    
def compare_metaprograms():
    omics_models_meta = load_file('model_metadata.csv', index_col=0)
    all_mmp_scores = load_file('celligner_geneset_scores.csv', local_dir=PROCESSED_DIR, index_col=0)
    geneset_table = load_file('geneset_table_updated.csv')
    
    # Gavish
    metaprograms = geneset_table[geneset_table['Source'] == "Gavish"]['Geneset'].unique().tolist()
    
    # run t-tests and Mann-Whitney U tests between 3D and 2D models in DepMap by lineage
    organoid_lineages = omics_models_meta[omics_models_meta['GrowthPattern'].isin(['Dome'])].value_counts('OncotreeLineage').index.tolist()
    organoid_lineages  = list(set(organoid_lineages) - {"Other", "Head and Neck", "Lung"})
    depmap_mmp_scores = all_mmp_scores.reindex(index=omics_models_meta.index.tolist(), columns=metaprograms).dropna(how='all', axis=0)
    
    # run t-tests and Mann-Whitney U tests between 3D and 2D models in DepMap by lineage
    organoid_vs_2d_mmp_lineage_tests = pd.concat([
        run_mwu_and_ttest(
            depmap_mmp_scores,
            omics_models_meta[
                (omics_models_meta['IsNextGen']) & 
                (omics_models_meta['OncotreeLineage'] == lin)
            ].index.tolist(),
            omics_models_meta[
                (omics_models_meta['IsTraditional2D']) & 
                (omics_models_meta['OncotreeLineage'] == lin)
            ].index.tolist(),
            group1_name='NextGen', group2_name='Traditional', nan_policy='omit'
        # ).assign(OncotreeLineage=lin) for lin in next_gen_lineages
        ).assign(OncotreeLineage=lin) for lin in organoid_lineages  
    ], axis=0).reset_index().rename({'index': 'Metaprogram'}, axis=1)
    
    # run t-tests and Mann-Whitney U tests between NextGen and 2D CNS/Brain models in DepMap
    neurosphere_vs_2d_mmp_lineage_tests = run_mwu_and_ttest(
        depmap_mmp_scores,
        omics_models_meta[
            (omics_models_meta['IsNextGen']) & 
            (omics_models_meta['OncotreeLineage'] == 'CNS/Brain')
        ].index.tolist(),
        omics_models_meta[
            (omics_models_meta['IsTraditional2D']) & 
            (omics_models_meta['OncotreeLineage'] == 'CNS/Brain')
        ].index.tolist(),
        group1_name='NextGen', group2_name='Traditional', nan_policy='omit'
    ).assign(OncotreeLineage='CNS/Brain').reset_index().rename({'index': 'Metaprogram'}, axis=1)
    
    # repeat for all Organoid lineages vs the 2D models from the same lineages
    organoid_vs_2d_mmp_test = run_mwu_and_ttest(
        depmap_mmp_scores,
        omics_models_meta[
            (omics_models_meta['IsNextGen']) &
            (omics_models_meta['OncotreeLineage'].isin(organoid_lineages))
        ].index.tolist(),
        omics_models_meta[
            (omics_models_meta['IsTraditional2D']) & 
            (omics_models_meta['OncotreeLineage'].isin(organoid_lineages))
        ].index.tolist(),
        group1_name='NextGen', group2_name='Traditional', nan_policy='omit'
    ).assign(OncotreeLineage='All Organoid Lineages').reset_index().rename({'index': 'Metaprogram'}, axis=1)
    
    # combine all the tests and reformat the missing data
    nextgen_vs_2d_mmp_tests = pd.concat([
        organoid_vs_2d_mmp_lineage_tests,
        neurosphere_vs_2d_mmp_lineage_tests,
        organoid_vs_2d_mmp_test
    ], axis=0, ignore_index=True)
    nextgen_vs_2d_mmp_tests['t p'] = nextgen_vs_2d_mmp_tests['t p'].fillna(1)
    nextgen_vs_2d_mmp_tests['t FDR'] = nextgen_vs_2d_mmp_tests['t FDR'].fillna(1)
    nextgen_vs_2d_mmp_tests['MWU p'] = nextgen_vs_2d_mmp_tests['MWU p'].fillna(1)
    nextgen_vs_2d_mmp_tests['MWU FDR'] = nextgen_vs_2d_mmp_tests['MWU FDR'].fillna(1)
    
    nextgen_vs_2d_mmp_tests.rename({
        'NextGen Mean': 'NextGenMean', 'NextGen n': 'NextGenN', 'Traditional Mean': 'TraditionalMean', 'Traditional n': 'TraditionalN',
        'Mean Difference': 'MeanDifference', 'MWU p': 'MannWhitneyPValue', 'MWU FDR': 'MannWhitneyFDR'
    }, axis=1).loc[
        :, [
            'Metaprogram', 'OncotreeLineage', 
            'NextGenMean', 'NextGenN', 'TraditionalMean', 'TraditionalN', 
            'MeanDifference', 'MannWhitneyPValue', 'MannWhitneyFDR'
        ]
    ].to_csv(os.path.join(PROCESSED_DIR, 'next_gen_vs_2d_metaprogram_expression_tests.csv'), index=False)
        
### Compute celligner nearest-neighbor metrics ### 

def identify_nearest_neighbors(dist_matrix, input_table, target_table, nn=25):
    top_nns = dist_matrix.loc[
        target_table.index.tolist(),
        input_table.index.tolist()
    ].apply(lambda x: pd.Series(x.nsmallest(nn).index.tolist()))
    
    top_nns = top_nns.melt(ignore_index=False).reset_index().rename(
        {'index': 'rank', 'variable': 'input', 'value': 'target'}, axis=1
    )
    top_nns['rank'] = top_nns['rank'] + 1
    
    top_nns = top_nns.merge(
        input_table.add_prefix('input_'), left_on='input', right_index=True
    ).merge(
        target_table.add_prefix('target_'), left_on='target', right_index=True
    )
    
    for c in target_table:
        top_nns[f'n_input_{c}_in_target'] = target_table.value_counts(c).reindex(
            index=top_nns[f'input_{c}'].tolist()
        ).fillna(0).values
    
    return top_nns.sort_index()

def classify_by_nearest_neighbors(nearest_neighbor_table, target_column='lineage', 
                                  top_n=25, exclude_insufficient=True):
    if top_n is None:
        top_n = nearest_neighbor_table['rank'].max()
        
    nn_table_subset = nearest_neighbor_table[
        (nearest_neighbor_table['rank'] <= top_n)
    ]
    
    if exclude_insufficient:
        nn_table_subset = nn_table_subset[
            (nn_table_subset[f'n_input_{target_column}_in_target']) >= (top_n / 2)
        ]
    
    confusion_count_matrix = nn_table_subset[['input', f'target_{target_column}']].groupby('input').value_counts().rename(
        'count'
    ).reset_index().pivot(
        index='input', columns=f'target_{target_column}', values='count'
    ).fillna(0)
    
    majority_predicted_class = confusion_count_matrix.idxmax(axis=1)
    majority_predicted_class_count = confusion_count_matrix.max(axis=1)
    neighbors_considered = confusion_count_matrix.sum(axis=1)
    true_class = nearest_neighbor_table.drop_duplicates(subset=['input', f'input_{target_column}']).set_index('input')[f'input_{target_column}']
    n_true_class = nearest_neighbor_table.drop_duplicates(subset=['input', f'n_input_{target_column}_in_target']).set_index('input')[f'n_input_{target_column}_in_target']
    
    predictions = pd.concat([
        majority_predicted_class.rename('predicted_class'),
        majority_predicted_class_count.rename('predicted_class_count'),
        neighbors_considered.rename('neighbors_considered'),
        true_class.rename('true_class'),
        n_true_class.rename('n_true_class_in_dataset')
    ], axis=1)
    
    return predictions

def run_celligner_classification():
    # Celligner data - import and map the gene ids to the right format
    celligner_distance_matrix = load_file('celligner_distance_matrix_w_hcmi_filtered.csv', index_col=0)
    omics_models_meta = load_file('model_metadata.csv', index_col=0)
    all_celligner_samples = load_file('celligner_coordinates_w_hcmi_filtered.csv', index_col=0)
    
    # partition the celligner table by dataset and model type

    # subset to the tcga table
    all_tcga_tumors_celligner_table = all_celligner_samples[
        all_celligner_samples['type'] == "TCGA+ tumor"
    ]

    # subset to the met500 table
    all_met500_tumors_celligner_table = all_celligner_samples[
        all_celligner_samples['type'] == 'MET500 tumor'
    ]

    # combine to an all-tumor table
    all_tumors_celligner_table = pd.concat([
        all_tcga_tumors_celligner_table,
        all_met500_tumors_celligner_table
    ], axis=0)

    # subset to the depmap models available for use first
    all_depmap_models_celligner_table = all_celligner_samples[
        all_celligner_samples.ModelID.isin(
            omics_models_meta[omics_models_meta["HasExprData"]].index.tolist())
    ]

    condition_to_model_map = all_depmap_models_celligner_table["ModelID"].to_dict()
    celligner_distance_matrix.rename(
        columns = condition_to_model_map, 
        inplace = True
    )

    all_depmap_models_celligner_table.set_index("ModelID", inplace = True)

    celligner_top_nns = identify_nearest_neighbors(
        celligner_distance_matrix, 
        input_table=all_depmap_models_celligner_table[['lineage', 'GrowthPattern']],
        target_table=all_tumors_celligner_table[['lineage']],
        nn=100
    )

    nextgen_dict = omics_models_meta["IsNextGen"].to_dict()
    traditional_dict = omics_models_meta["IsTraditional2D"].to_dict()
    adherent_dict = omics_models_meta["IsAdherent2D"].to_dict()
    
    celligner_top_nns['inputIsNextGen'] = celligner_top_nns["input"].map(nextgen_dict)
    celligner_top_nns['inputIsTraditional2D'] = celligner_top_nns["input"].map(traditional_dict)
    celligner_top_nns['inputIsAdherent2D'] = celligner_top_nns["input"].map(adherent_dict)

    celligner_predictions_lineage = classify_by_nearest_neighbors(
        celligner_top_nns, top_n=25, target_column='lineage', exclude_insufficient=False
    ).join(omics_models_meta[['GrowthPattern', 'IsNextGen', 'IsTraditional2D', 'IsAdherent2D']]).reset_index().rename(
        {'input': 'ModelID'}, axis=1
    )
    celligner_predictions_lineage['immediate_fail'] = celligner_predictions_lineage['n_true_class_in_dataset'] < (25 / 2)

    # perform lineage prediction for novartix pdx subset
    all_pdx_models_celligner_table = all_celligner_samples[
        all_celligner_samples["type"] == "Novartis_PDX"
    ]
    celligner_top_nns_pdx = identify_nearest_neighbors(
        celligner_distance_matrix, 
        input_table= all_pdx_models_celligner_table[["lineage"]], 
        target_table= all_tumors_celligner_table[["lineage"]]
    )
    celligner_predictions_lineage_pdx = classify_by_nearest_neighbors(
        celligner_top_nns_pdx, top_n= 25, target_column= "lineage", exclude_insufficient= False
    ).reset_index().rename({'input': 'ModelID'}, axis=1)
    celligner_predictions_lineage_pdx["immediate_fail"] = celligner_predictions_lineage_pdx["n_true_class_in_dataset"] < (25 / 2)
    
    celligner_top_nns.to_csv(os.path.join(PROCESSED_DIR, 'celligner_nearest_neighbor_table.csv'), index=False)
    celligner_predictions_lineage.to_csv(os.path.join(PROCESSED_DIR, 'celligner_lineage_predictions.csv'), index=False)
    celligner_top_nns_pdx.to_csv(os.path.join(PROCESSED_DIR, 'celligner_nearest_neighbor_table_pdx.csv'), index=False)
    celligner_predictions_lineage_pdx.to_csv(os.path.join(PROCESSED_DIR, 'celligner_lineage_predictions_pdx.csv'), index=False)
    
##### Dependency Biomarker Analysis #####

# Expression addictions

def compute_expression_addictions(screen_metadata_subset, screen_gene_effect, screen_expr):
    screen_effects_matrix, screen_expr_matrix = screen_gene_effect.reindex(
        index=screen_metadata_subset.index.tolist()
    ).dropna(how='all', axis=1).align(screen_expr, join='inner')
    self_expression_tests = run_correlations_paired(
        [(x, x) for x in screen_effects_matrix.columns.tolist()], screen_effects_matrix, screen_expr_matrix
    ).rename({'Feature 1': 'Dependency', 'Feature 2': 'Biomarker'}, axis=1)
    return self_expression_tests

def run_expression_addiction_analysis():
    screen_gene_effect = load_file('screen_gene_effect.csv', index_col=0)
    screen_metadata = load_file('screen_metadata.csv', index_col=0)
    expr = load_full_matrix('expression')
    screen_expr = expand_model_matrix_to_screens(expr, screen_metadata.reset_index().set_index('ModelID')['ScreenID'])
    
    next_gen_screens = screen_metadata[
        screen_metadata['ScreenType'].isin(['2DO', '3DO', '3DN', '2DN']) &
        screen_metadata['IsNextGen'] & 
        (screen_metadata['OncotreeLineage'] != 'Other')
    ]
    next_gen_lineages = screen_metadata[screen_metadata['IsNextGen']].value_counts('OncotreeLineage').index.tolist()
    
    adherent_next_gen_lineage_matched_screens = screen_metadata[
        screen_metadata['IsTraditional2D'] & 
        screen_metadata['OncotreeLineage'].isin(next_gen_lineages)
    ]
    
    next_gen_self_expression_tests = compute_expression_addictions(next_gen_screens, screen_gene_effect, screen_expr)
    adherent_lineage_matched_self_expression_tests = compute_expression_addictions(adherent_next_gen_lineage_matched_screens, screen_gene_effect, screen_expr)
    
    expression_addiction_tests = next_gen_self_expression_tests.merge(
        adherent_lineage_matched_self_expression_tests.loc[:, ['Biomarker', 'Pearson Correlation']].rename(
            {'Pearson Correlation': 'TraditionalPearsonR'}, axis=1), 
        left_on=['Biomarker'], right_on=['Biomarker']
    )
    expression_addiction_tests = expression_addiction_tests.rename({
        'Pearson Correlation': 'NextGenPearsonR', 'p-value': 'NextGenPValue', 'n': 'NextGenN', 
    }, axis=1)
    
    expression_addiction_tests.to_csv(os.path.join(PROCESSED_DIR, 'expression_addiction_tests_full.csv'), index=False)

# Paralog dependencies

def run_paralog_dependency_analysis(paralog_table, gene_effect_matrix, biomarker_matrix, biomarker_type):
    paralog_input_table = paralog_table.copy()
    
    paralog_input_table['Gene CRISPR Var'] = gene_effect_matrix.var().reindex(
        index=paralog_input_table['Gene CDS ID']
    ).values
    paralog_input_table[f'Paralog {biomarker_type} Var'] = biomarker_matrix.var().reindex(
        index=paralog_input_table['Paralog CDS ID']
    ).values
    
    # ignore anything without variation in the dependency or biomarker
    dependencies_to_paralogs = paralog_input_table[
        (paralog_input_table['Gene CRISPR Var'] > 0) &
        (paralog_input_table[f'Paralog {biomarker_type} Var'] > 0)
    ].groupby('Gene CDS ID').apply(lambda x: list(x['Paralog CDS ID'])).to_dict()
    
    # for each dependency, run the analysis on all its paralogs
    paralog_associations = []
    for dep in tqdm(dependencies_to_paralogs):
        expected_biomarkers = dependencies_to_paralogs[dep]
        biomarker_tests = get_preselected_biomarkers(
            dep, 
            expected_biomarkers,
            sample_subset=gene_effect_matrix.index.tolist(),
            gene_effect_matrix=gene_effect_matrix,
            biomarker_matrix=biomarker_matrix
        )
        biomarker_tests['Dependency'] = dep
        paralog_associations.append(biomarker_tests)
    paralog_associations = pd.concat(paralog_associations, axis=0, ignore_index=True)
    paralog_associations['Biomarker Type'] = biomarker_type
    
    paralog_associations[f'Biomarker Var'] = biomarker_matrix.var().reindex(
        index=paralog_associations['Biomarker']
    ).values
    
    return paralog_associations.loc[
        :, ['Dependency', 'Biomarker', 'Biomarker Type', 'Biomarker Var', 'Pearson Correlation', 'p-value', 'n']
    ]

def compute_paralog_dependencies(screen_metadata_subset, ensembl_paralogs, screen_gene_effect, screen_expr):
    screen_effects_matrix = screen_gene_effect.reindex(index=screen_metadata_subset.index.tolist())
    screen_expr_matrix = screen_expr.reindex(index=screen_effects_matrix.index.tolist())
    
    paralog_associations = run_paralog_dependency_analysis(
        ensembl_paralogs, screen_effects_matrix, screen_expr_matrix, biomarker_type='EXPR'
    )
    
    return paralog_associations

def run_paralog_analysis():
    screen_gene_effect = load_file('screen_gene_effect.csv', index_col=0)
    screen_metadata = load_file('screen_metadata.csv', index_col=0)
    expr = load_full_matrix('expression')
    screen_expr = expand_model_matrix_to_screens(expr, screen_metadata.reset_index().set_index('ModelID')['ScreenID'])
    ensembl_paralogs = load_file('paralog_pair_table.csv')
    
    next_gen_screens = screen_metadata[
        screen_metadata['ScreenType'].isin(['2DO', '3DO', '3DN', '2DN']) &
        screen_metadata['IsNextGen'] & 
        (screen_metadata['OncotreeLineage'] != 'Other')
    ]
    next_gen_lineages = screen_metadata[screen_metadata['IsNextGen']].value_counts('OncotreeLineage').index.tolist()
    
    adherent_next_gen_lineage_matched_screens = screen_metadata[
        screen_metadata['IsTraditional2D'] & 
        screen_metadata['OncotreeLineage'].isin(next_gen_lineages)
    ]
    
    next_gen_paralog_associations = compute_paralog_dependencies(
        next_gen_screens, ensembl_paralogs, screen_gene_effect, screen_expr
    )
    adherent_lin_matched_paralog_associations = compute_paralog_dependencies(
        adherent_next_gen_lineage_matched_screens, ensembl_paralogs, screen_gene_effect, screen_expr
    )
    
    paralog_dependency_tests = next_gen_paralog_associations.merge(
        adherent_lin_matched_paralog_associations.loc[:, ['Biomarker', 'Dependency', 'Pearson Correlation']].rename(
            {'Pearson Correlation': 'TraditionalPearsonR'}, axis=1), 
        left_on=['Biomarker', 'Dependency'], right_on=['Biomarker', 'Dependency'], how='left'
    )
    paralog_dependency_tests = paralog_dependency_tests.rename({
        'Biomarker Type': 'BiomarkerType', 'Biomarker Var': 'BiomarkerVariance', 
        'Pearson Correlation': 'NextGenPearsonR', 'p-value': 'NextGenPValue', 'n': 'NextGenN'
    }, axis=1)
    
    paralog_dependency_tests.to_csv(os.path.join(PROCESSED_DIR, 'paralog_dependency_tests_full.csv'), index=False)

# Oncogene / tumor suppressor alteration-associated dependencies
    
def prepare_alteration_data(amplification_min_cn = 3, deletion_max_cn = 0.25):
    # load gene annotations
    geneset_table = load_file('geneset_table_updated.csv')
    og_genes = geneset_table[(geneset_table['Geneset'] == "Oncogenes")]['cds_gene_id'].tolist()
    tsg_genes = geneset_table[(geneset_table['Geneset'] == "Tumor suppressors")]['cds_gene_id'].tolist()
    
    # load crispr and omics data
    screen_gene_effect = load_file('screen_gene_effect.csv', index_col=0)
    screen_gene_dependency = load_file('screen_gene_dependency.csv', index_col=0)
    expr = load_full_matrix('expression')
    cn = load_full_matrix('copy_number')
    hotspot = load_full_matrix('hotspot')
    damaging = load_full_matrix('damaging')
    
    # load screen sets
    screen_metadata = load_file('screen_metadata.csv', index_col=0)
    next_gen_screens = screen_metadata[
        screen_metadata['ScreenType'].isin(['2DO', '3DO', '3DN', '2DN']) &
        screen_metadata['IsNextGen'] &
        (screen_metadata['OncotreeLineage'] != 'Other')
    ]
    next_gen_screens = next_gen_screens.loc[
        next_gen_screens.index.isin(screen_gene_effect.index)
    ]
    next_gen_lineages = screen_metadata[screen_metadata['IsNextGen']].value_counts('OncotreeLineage').index.tolist()
    adherent_next_gen_lineage_matched_screens = screen_metadata[
        screen_metadata['IsTraditional2D'] & 
        screen_metadata['OncotreeLineage'].isin(next_gen_lineages)
    ]
    adherent_next_gen_lineage_matched_screens = adherent_next_gen_lineage_matched_screens.loc[
        adherent_next_gen_lineage_matched_screens.index.isin(screen_gene_effect.index)
    ]
    
    # reformat omics data to match crispr
    screen_expr = expand_model_matrix_to_screens(expr, screen_metadata.reset_index().set_index('ModelID')['ScreenID'])
    screen_cn = expand_model_matrix_to_screens(cn, screen_metadata.reset_index().set_index('ModelID')['ScreenID'])
    screen_hotspot = expand_model_matrix_to_screens(hotspot, screen_metadata.reset_index().set_index('ModelID')['ScreenID'])
    screen_damaging = expand_model_matrix_to_screens(damaging, screen_metadata.reset_index().set_index('ModelID')['ScreenID'])
    
    # identify all onc/tsgs with at least one hotspot mutation in the next gen models
    either_with_hotspot = pd.DataFrame(((screen_hotspot > 0).reindex(
        index=next_gen_screens.index.tolist(),
        columns=og_genes + tsg_genes
    ) > 0).sum().loc[lambda x: x > 0].rename('NextGenHotspotModels'))
    either_with_hotspot['Oncogene'] = either_with_hotspot.index.isin(og_genes)
    either_with_hotspot['TumorSuppressor'] = either_with_hotspot.index.isin(tsg_genes)
    either_with_hotspot['GeneSymbol'] = either_with_hotspot.index.str.split(' ').str[0]
    
    cancer_gene_annotations = pd.concat([
        pd.Series(True, index=og_genes).rename('Oncogene'),
        pd.Series(True, index=tsg_genes).rename('TumorSuppressor')
    ], axis=1).fillna(False)
    
    # combine matrices to mark the presence of oncogenic or tummor suppressive alterations
    nextgen_hotspot_matrix = screen_hotspot.reindex(
        index=next_gen_screens.index.tolist(),
        columns=cancer_gene_annotations.index.tolist()
    ).dropna(how='all', axis=0) > 0
    nextgen_damaging_matrix = screen_damaging.reindex(
        index=next_gen_screens.index.tolist(),
        columns=cancer_gene_annotations.index.tolist()
    ).dropna(how='all', axis=0) > 0
    nextgen_cn_matrix = screen_cn.reindex(
        index=next_gen_screens.index.tolist(),
        columns=cancer_gene_annotations.index.tolist()
    ).dropna(how='all', axis=0)
    nextgen_amplified_mtx = nextgen_cn_matrix > amplification_min_cn
    nextgen_deleted_mtx = nextgen_cn_matrix < deletion_max_cn
    nextgen_og_altered_mtx = nextgen_hotspot_matrix | nextgen_amplified_mtx
    nextgen_ts_altered_mtx = nextgen_damaging_matrix | nextgen_deleted_mtx
    
    # summarize the occurrence of alterations in the NextGen dataset
    cancer_gene_annotations['AnyHotspotInNextGen'] = cancer_gene_annotations.index.isin(either_with_hotspot.index.tolist())
    cancer_gene_annotations['n_hotspots'] = nextgen_hotspot_matrix.sum().reindex(cancer_gene_annotations.index.tolist()).fillna(0)
    cancer_gene_annotations['n_damaging'] = nextgen_damaging_matrix.sum().reindex(cancer_gene_annotations.index.tolist()).fillna(0)
    cancer_gene_annotations['n_amplified'] = nextgen_amplified_mtx.sum().reindex(cancer_gene_annotations.index.tolist()).fillna(0)
    cancer_gene_annotations['n_deleted'] = nextgen_deleted_mtx.sum().reindex(cancer_gene_annotations.index.tolist()).fillna(0)
    cancer_gene_annotations['n_oncogenic_alteration'] = nextgen_og_altered_mtx.sum().reindex(cancer_gene_annotations.index.tolist()).fillna(0)
    cancer_gene_annotations['n_tumor_suppressor_alteration'] = nextgen_ts_altered_mtx.sum().reindex(cancer_gene_annotations.index.tolist()).fillna(0)
    
    # only consider dependency genes which are vulnerabilities in at least 3 models
    nextgen_is_screen_dependency = screen_gene_dependency.loc[next_gen_screens.index.tolist(), :] > 0.5
    nextgen_is_screen_dependency = nextgen_is_screen_dependency.loc[
        :, nextgen_is_screen_dependency.sum().loc[lambda x: x >= 3].index.tolist()
    ]
    next_gen_screen_effects = screen_gene_effect.reindex_like(nextgen_is_screen_dependency)
    
    # repeat for adherent
    # combine matrices to mark the presence of oncogenic or tummor suppressive alterations
    adherent_hotspot_matrix = screen_hotspot.reindex(
        index=adherent_next_gen_lineage_matched_screens.index.tolist(),
        columns=cancer_gene_annotations.index.tolist()
    ).dropna(how='all', axis=0) > 0
    adherent_damaging_matrix = screen_damaging.reindex(
        index=adherent_next_gen_lineage_matched_screens.index.tolist(),
        columns=cancer_gene_annotations.index.tolist()
    ).dropna(how='all', axis=0) > 0
    adherent_cn_matrix = screen_cn.reindex(
        index=adherent_next_gen_lineage_matched_screens.index.tolist(),
        columns=cancer_gene_annotations.index.tolist()
    ).dropna(how='all', axis=0)
    adherent_amplified_mtx = adherent_cn_matrix > 3
    adherent_deleted_mtx = adherent_cn_matrix < 0.25
    adherent_og_altered_mtx = adherent_hotspot_matrix | adherent_amplified_mtx
    adherent_ts_altered_mtx = adherent_damaging_matrix | adherent_deleted_mtx
    
    onc_gof_altered_mtx = pd.concat([
        nextgen_og_altered_mtx,
        adherent_og_altered_mtx
    ], axis=0)
    tsg_lof_altered_mtx = pd.concat([
        nextgen_ts_altered_mtx,
        adherent_ts_altered_mtx
    ], axis=0)
    
    onc_gof_altered_mtx.astype(float).to_csv(os.path.join(PROCESSED_DIR, 'oncogene_gof_alteration_matrix.csv'), index=True)
    tsg_lof_altered_mtx.astype(float).to_csv(os.path.join(PROCESSED_DIR, 'tumor_suppressor_lof_alteration_matrix.csv'), index=True)
    
    # match the adherent features to the same features tested in nextgen
    adherent_is_screen_dependency = (
        screen_gene_dependency.loc[adherent_next_gen_lineage_matched_screens.index.tolist(), :] > 0.5
    ).reindex(columns=nextgen_is_screen_dependency.columns.tolist())
    adherent_screen_effects = screen_gene_effect.reindex_like(adherent_is_screen_dependency)
    
    return (
        cancer_gene_annotations, screen_gene_effect, screen_metadata,
        next_gen_screen_effects, nextgen_og_altered_mtx, nextgen_ts_altered_mtx,
        adherent_screen_effects, adherent_og_altered_mtx, adherent_ts_altered_mtx
    )

# run t-tests and Mann-Whitney U tests
def run_oncogene_tsg_dependency_analysis(genes_to_test, screen_effect_matrix, alteration_matrix, dependency_matrix=None):
    association_tests = []
    for g in tqdm(genes_to_test):
        biomarker_tests = run_mwu_and_ttest(
            screen_effect_matrix,
            group1=alteration_matrix[g].loc[lambda x: x].index.tolist(),
            group1_name="Altered",
            group2=alteration_matrix[g].loc[lambda x: ~x].index.tolist(),
            group2_name='WT',
            fishers_exact_matrix=dependency_matrix,
            nan_policy='omit'
        ).reset_index().rename({'index': 'Dependency'}, axis=1).assign(Biomarker=g).drop(
            ['t FDR', 'U1', 'MWU FDR', 'Fisher FDR'], axis=1, errors='ignore'
        )
        association_tests.append(biomarker_tests)
        
    association_tests = pd.concat(association_tests, axis=0, ignore_index=True)
    association_tests = association_tests.loc[
        :, ['Dependency', 'Biomarker'] + [
            x for x in association_tests.columns.tolist() if x not in ['Dependency', 'Biomarker']
        ]
    ]
    return association_tests

def run_dependency_vs_lineage_biomarker_analysis(screen_effect_matrix, categorical_annotation, dependency_matrix=None):
    category_to_associations = dict()
    categorical_annotation = categorical_annotation.loc[screen_effect_matrix.index.tolist()].copy()
    categories = categorical_annotation.unique().tolist()
    for c in tqdm(categories):
        category_to_associations[c] = run_mwu_and_ttest(
            screen_effect_matrix,
            group1=categorical_annotation.loc[lambda x: x == c].index.tolist(),
            group1_name=f"Ingroup",
            group2=categorical_annotation.loc[lambda x: x != c].index.tolist(),
            group2_name=f"Outgroup",
            fishers_exact_matrix=dependency_matrix,
            nan_policy='omit'
        ).reset_index().rename({'index': 'Dependency'}, axis=1)
    
    return pd.concat(category_to_associations, axis=0).reset_index().rename({'level_0': 'OncotreeLineage'}, axis=1).drop('level_1', axis=1)

def run_oncogene_tsg_association_analysis():
    (
        cancer_gene_annotations, screen_gene_effect, screen_metadata,
        next_gen_screen_effects, nextgen_og_altered_mtx, nextgen_ts_altered_mtx,
        adherent_screen_effects, adherent_og_altered_mtx, adherent_ts_altered_mtx
    ) = prepare_alteration_data()
    
    # select the genes to test for each gene type
    oncogenes_to_test = cancer_gene_annotations[
        (cancer_gene_annotations['Oncogene']) & (cancer_gene_annotations['n_oncogenic_alteration'] >= 5)
    ].index.tolist()
    tsgs_to_test = cancer_gene_annotations[
        (cancer_gene_annotations['TumorSuppressor']) & (cancer_gene_annotations['n_tumor_suppressor_alteration'] >= 5)
    ].index.tolist()

    # run analysis on oncogenes in next gen models
    next_gen_oncogene_associations = run_oncogene_tsg_dependency_analysis(
        oncogenes_to_test, next_gen_screen_effects, nextgen_og_altered_mtx
    )
    next_gen_oncogene_associations = add_cytobands(next_gen_oncogene_associations, column_names=['Dependency', 'Biomarker'])
    next_gen_oncogene_associations['AlterationType'] = 'Oncogene'

    # run analysis on tumor suppressors in next gen models
    next_gen_tsg_associations = run_oncogene_tsg_dependency_analysis(
        tsgs_to_test, next_gen_screen_effects, nextgen_ts_altered_mtx
    )
    next_gen_tsg_associations = add_cytobands(next_gen_tsg_associations, column_names=['Dependency', 'Biomarker'])
    next_gen_tsg_associations['AlterationType'] = 'Tumor Suppressor'
    
    # merge the gene types
    next_gen_oncogene_associations = next_gen_oncogene_associations.merge(cancer_gene_annotations['AnyHotspotInNextGen'], left_on = 'Biomarker', right_index=True)
    next_gen_tsg_associations = next_gen_tsg_associations.merge(cancer_gene_annotations['AnyHotspotInNextGen'], left_on = 'Biomarker', right_index=True)

    next_gen_onc_tsg_associations = pd.concat([
        next_gen_oncogene_associations,
        next_gen_tsg_associations
    ], axis=0, ignore_index=True)
        
    # run analysis on oncogenes in adherent models
    adherent_oncogene_associations = run_oncogene_tsg_dependency_analysis(
        oncogenes_to_test, adherent_screen_effects, adherent_og_altered_mtx
    )
    adherent_oncogene_associations = add_cytobands(adherent_oncogene_associations, column_names=['Dependency', 'Biomarker'])
    adherent_oncogene_associations['AlterationType'] = 'Oncogene'

    # run analysis on tumor suppressors in adherent models
    adherent_tsg_associations = run_oncogene_tsg_dependency_analysis(
        tsgs_to_test, adherent_screen_effects, adherent_ts_altered_mtx
    )
    adherent_tsg_associations = add_cytobands(adherent_tsg_associations, column_names=['Dependency', 'Biomarker'])
    adherent_tsg_associations['AlterationType'] = 'Tumor Suppressor'
    
    # merge the gene types
    adherent_oncogene_associations = adherent_oncogene_associations.merge(cancer_gene_annotations['AnyHotspotInNextGen'], left_on = 'Biomarker', right_index=True)
    adherent_tsg_associations = adherent_tsg_associations.merge(cancer_gene_annotations['AnyHotspotInNextGen'], left_on = 'Biomarker', right_index=True)

    adherent_onc_tsg_associations = pd.concat([
        adherent_oncogene_associations,
        adherent_tsg_associations
    ], axis=0, ignore_index=True)
    
    onc_tsg_associations_full = next_gen_onc_tsg_associations.merge(
        adherent_onc_tsg_associations.loc[:, ['Dependency', 'Biomarker', 'Mean Difference', 'MWU p']].rename({
            'Mean Difference': 'TraditionalMeanDifference', 'MWU p': 'TraditionalMannWhitneyPValue'
        }, axis=1), how='left', left_on=['Dependency', 'Biomarker'], right_on=['Dependency', 'Biomarker'])
    onc_tsg_associations_full = onc_tsg_associations_full.rename({
        'Altered Mean': 'NextGenAlteredMean', 'Altered Var': 'NextGenAlteredVariance', 'Altered n': 'NextGenAlteredN', 
        'WT Mean': 'NextGenWildTypeMean', 'WT Var': 'NextGenWildTypeVariance', 'WT n': 'NextGenWildTypeN', 
        'Mean Difference': 'NextGenMeanDifference', 'MWU p': 'NextGenMannWhitneyPValue', 'AlterationType': 'BiomarkerAlterationType'
    }, axis=1).loc[
        :, [
            'Dependency', 'Biomarker', 'BiomarkerAlterationType',
            'NextGenAlteredMean', 'NextGenAlteredVariance', 'NextGenAlteredN', 
            'NextGenWildTypeMean', 'NextGenWildTypeVariance', 'NextGenWildTypeN',
            'NextGenMeanDifference', 'NextGenMannWhitneyPValue', 
            'TraditionalMeanDifference', 'TraditionalMannWhitneyPValue'
        ]
    ]
    onc_tsg_associations_full.to_csv(os.path.join(PROCESSED_DIR, 'onc_tsg_alteration_dependency_tests_full.csv'), index=False)
    
    # add tests for lineage
    nextgen_dependency_lineage_tests = run_dependency_vs_lineage_biomarker_analysis(
        next_gen_screen_effects, 
        screen_metadata['OncotreeLineage'], 
    )
    nextgen_dependency_lineage_tests.to_csv(os.path.join(PROCESSED_DIR, 'next_gen_dependency_lineage_tests.csv'), index=False)
    
    cancer_gene_annotations.rename_axis('Gene').to_csv(os.path.join(PROCESSED_DIR, 'next_gen_cancer_gene_annotations.csv'))
    
# Glioblastoma dependency analysis

def compute_gbm_cdk6_dependency_biomarker_tests(correlation_threshold=0.3, correlation_significance_threshold=0.005,
                                                 min_mutations=5, mean_effect_threshold=0.3, mwu_significance_threshold=0.05):
    screen_metadata = load_file('screen_metadata.csv', index_col=0)
    screen_gene_effect = load_file('screen_gene_effect.csv', index_col=0)
    expr = load_full_matrix('expression')
    hotspot = load_full_matrix('hotspot')
    damaging = load_full_matrix('damaging')

    gbm_screens_expr = screen_metadata[
        screen_metadata['OncotreeSubtype'].fillna('NA').str.contains('Glioblastoma|Gliosarcoma') &
        screen_metadata['HasExprData']
    ].index.tolist()
    
    # reformat model-indexed data into screen-indexed
    screen_expr = expand_model_matrix_to_screens(expr, screen_metadata.reset_index().set_index('ModelID')['ScreenID'])
    screen_hotspot = expand_model_matrix_to_screens(hotspot, screen_metadata.reset_index().set_index('ModelID')['ScreenID'])
    screen_damaging = expand_model_matrix_to_screens(damaging, screen_metadata.reset_index().set_index('ModelID')['ScreenID'])
    screen_mutation_matrix = ((screen_hotspot > 0) | (screen_damaging.reindex_like(hotspot) > 0))
    # ID expression-based biomarkers
    expr_variance = screen_expr.reindex(index=gbm_screens_expr).var()
    cdk6_dep_expr_correlations = quick_corr(
        screen_gene_effect.reindex(index=gbm_screens_expr)[search_gene('CDK6')],
        screen_expr.reindex(index=gbm_screens_expr).loc[:, expr_variance > 0], n=None
    )
    cdk6_dep_expr_correlations['ExprVariance'] = expr_variance.reindex(index=cdk6_dep_expr_correlations.index.tolist())
    cdk6_expr_volcano_df = cdk6_dep_expr_correlations[cdk6_dep_expr_correlations['ExprVariance'] > 1].copy()
    cdk6_expr_volcano_df['FDR'] = scipy.stats.false_discovery_control(cdk6_expr_volcano_df['p-value'])
    cdk6_expr_volcano_df['Significant'] = (
        (cdk6_expr_volcano_df['Pearson Correlation'].abs() > correlation_threshold) &
        (cdk6_expr_volcano_df['FDR'] < correlation_significance_threshold)
    )
    cdk6_expr_volcano_df.rename_axis('Biomarker').to_csv(
        os.path.join(PROCESSED_DIR, 'gbm_cdk6_dependency_expression_tests.csv')
    )

    gbm_screens_cn = screen_metadata[
        screen_metadata['OncotreeSubtype'].fillna('NA').str.contains('Glioblastoma|Gliosarcoma')
        & screen_metadata['HasCNData']
    ].index.tolist()
    
    gbm_mutation_biomarkers = screen_mutation_matrix.reindex(index=gbm_screens_cn).sum().loc[lambda x: x >= min_mutations].index.tolist()
    gbm_cdk6_dep_mut_tests = []
    for bm in gbm_mutation_biomarkers:
        bm_series = screen_mutation_matrix.reindex(index=gbm_screens_cn)[bm].dropna()
        altered = bm_series[bm_series].index.tolist()
        wt = bm_series[~bm_series].index.tolist()
        cdk6_mut_test = run_mwu_and_ttest(
            screen_gene_effect.reindex(index=gbm_screens_cn)[[search_gene('CDK6')]],
            altered,
            wt,
            group1_name='Altered',
            group2_name='WT',
            nan_policy='omit'
        ).assign(Biomarker=bm)

        gbm_cdk6_dep_mut_tests.append(cdk6_mut_test)
    gbm_cdk6_dep_mut_tests = pd.concat(gbm_cdk6_dep_mut_tests, axis=0).rename_axis('Dependency').reset_index().loc[
        :, ['Dependency', 'Biomarker', 'Altered Mean', 'Altered n', 'WT Mean', 'WT n', 'Mean Difference', 'MWU p']
    ]
    gbm_cdk6_dep_mut_tests['FDR'] = scipy.stats.false_discovery_control(gbm_cdk6_dep_mut_tests['MWU p'])
    gbm_cdk6_dep_mut_tests['Significant'] = (
        (gbm_cdk6_dep_mut_tests['Mean Difference'].abs() > 0.3) &
        (gbm_cdk6_dep_mut_tests['FDR'] < 0.05)
    )
    gbm_cdk6_dep_mut_tests.to_csv(os.path.join(PROCESSED_DIR, 'gbm_cdk6_dependency_mutation_tests.csv'), index=False)
    
# PDAC-classical associated dependencies
    
def run_pdac_classical_analysis():
    # identify the models that will be used for analysis
    all_mmp_scores = load_file('celligner_geneset_scores.csv', local_dir=PROCESSED_DIR, index_col=0)
    screen_metadata = load_file('screen_metadata.csv', index_col=0)
    screen_gene_effect = load_file('screen_gene_effect.csv', index_col=0)
    all_mmp_scores.index = all_mmp_scores.index.astype(str)
    screen_gene_effect.index = screen_gene_effect.index.astype(str)

    organoid_screens = screen_metadata[
        screen_metadata["PassesQC"] &
        screen_metadata['IsNextGen'] &
        screen_metadata['ScreenType'].isin(['2DO', '3DO']) & 
        (screen_metadata['OncotreeLineage'] != 'Other') & 
        (screen_metadata['OncotreeLineage'] != 'CNS/Brain')
        & screen_metadata.index.isin(screen_gene_effect.index.tolist())
    ]    
    
    
    x = (all_mmp_scores.reindex(index=organoid_screens['ModelID'].tolist()).rename(screen_metadata.reset_index().set_index('ModelID')['ScreenID'])['PDAC-classical'])
    y = (screen_gene_effect.loc[organoid_screens.index.tolist(), :].dropna(how='all', axis=1))
    x.index = x.index.map(str)
    y.index = y.index.map(str)


    mmp_pdac_classical_associations = quick_corr(
        x, y, n=None
    )

    mmp_pdac_classical_associations = mmp_pdac_classical_associations.reset_index().rename(
        {'index': 'Dependency', 'Pearson Correlation': 'NextGenPearsonR', 'p-value': 'NextGenPValue', 'n': 'NextGenN'}, axis=1
    ).assign(Biomarker='PDAC-classical').loc[
        :, ['Dependency', 'Biomarker', 'NextGenPearsonR', 'NextGenPValue', 'NextGenN']
    ]
    mmp_pdac_classical_associations.to_csv(os.path.join(PROCESSED_DIR, 'pdac_classical_dependency_tests_full.csv'), index=False)
    
##### Global dependency differences across formats #####

def compare_next_gen_vs_adherent_dependency():
    screen_gene_effect = load_file('screen_gene_effect.csv', index_col=0)
    screen_gene_dependency = load_file('screen_gene_dependency.csv', index_col=0)
    screen_metadata = load_file('screen_metadata.csv', index_col=0)
    geneset_table = load_file('geneset_table_updated.csv')
    
    # identify the models that will be used for analysis
    organoid_screens = screen_metadata[
        screen_metadata['IsNextGen'] & 
        screen_metadata['ScreenType'].isin(["3DO"]) &
        (screen_metadata['OncotreeLineage'] != 'Other')
        & (screen_metadata.index.isin(screen_gene_effect.index.tolist()))
    ]
    
    adherent_organoid_lineage_matched_screens = screen_metadata[
        screen_metadata['IsAdherent2D'] &
        screen_metadata['OncotreeLineage'].isin(organoid_screens['OncotreeLineage'].unique().tolist())
        & screen_metadata.index.isin(screen_gene_effect.index.tolist())
    ]
    
    # calculate mean differences across organoid and adherent models
    global_organoids_vs_2d_difference_tests = run_mwu_and_ttest(
        screen_gene_effect,
        organoid_screens.index.tolist(),
        adherent_organoid_lineage_matched_screens.index.tolist(),
        group1_name='Organoids',
        group2_name='CellLinesOrganoidLineageMatched',
        nan_policy='omit'
    )
    
    # identify genes which should be excluded from analysis
    ess_screen_medians = (
        (screen_gene_effect.loc[organoid_screens.index.tolist(), :].median() < -1) & 
        (screen_gene_effect.loc[adherent_organoid_lineage_matched_screens.index.tolist(), :].median() < -1)
    ).rename('Screen GE Median < -1') # essential genes: strong depletion effects in either set

    ess_prob_dep = (
        ((screen_gene_dependency.loc[organoid_screens.index.tolist(), :] > 0.5).mean() > 0.5) &
        ((screen_gene_dependency.loc[adherent_organoid_lineage_matched_screens.index.tolist(), :] > 0.5).mean() > 0.5)
    ).rename('> 50% Prob of Dep in > 50% Lines') # essential genes: consistently dependent in either set

    # drop genes which are essential
    ess_to_remove = pd.concat([
        ess_screen_medians, 
        ess_prob_dep
    ], axis=1)
    ess_to_remove['drop_essential'] = (
        ess_to_remove['Screen GE Median < -1'] | 
        ess_to_remove['> 50% Prob of Dep in > 50% Lines']
    )
    
    # annotate genes with their effect size differences
    org_2d_diff_df = global_organoids_vs_2d_difference_tests.dropna(
        subset=['Organoids Mean', 'CellLinesOrganoidLineageMatched Mean']
    ).sort_values(
        'Mean Difference'
    ).join(ess_to_remove['drop_essential'])
    org_2d_diff_df = org_2d_diff_df.reset_index().rename({'index': 'Dependency'}, axis=1)
    org_2d_diff_df.rename({
        'Organoids Mean': 'OrganoidsMean', 'Organoids n': 'OrganoidsN', 'CellLinesOrganoidLineageMatched Mean': 'AdherentMean', 'CellLinesOrganoidLineageMatched n': 'AdherentN', 'Mean Difference': 'MeanDifference', 'MWU p': 'MannWhitneyPValue', 'MWU FDR': 'MannWhitneyFDR', 'drop_essential': 'DropEssential'
    }, axis=1).loc[
        :, [
            'Dependency', 'OrganoidsMean', 'OrganoidsN', 'AdherentMean', 'AdherentN', 'MeanDifference', 'MannWhitneyPValue', 'MannWhitneyFDR', 'DropEssential'
        ]
    ].to_csv(os.path.join(PROCESSED_DIR, 'organoids_vs_2d_gene_dependency_tests_full.csv'), index=False)
    
    # run gene set analysis on the distribution of dependency differences across model sets
    # load genesets
    hallmark = geneset_table[geneset_table['Source'] == "Hallmark"].groupby('Geneset').apply(lambda x: list(x['GeneSymbol'])).to_dict()
    kegg_legacy = geneset_table[geneset_table['Source'] == "KEGG"].groupby('Geneset').apply(lambda x: list(x['GeneSymbol'])).to_dict()

    # perform overrepresentation analysis in the tails
    rank_hypergeometric_results = dict()
    for collection, genesets in zip(['Hallmark', 'KEGG'], [hallmark, kegg_legacy]):
        rank_hypergeometric_results[collection] = run_rank_hypergeometric(
            org_2d_diff_df[~org_2d_diff_df['drop_essential']].set_index('Dependency')['Mean Difference'].sort_values(),
            genesets,
            n_top=200, direction='both',
            all_genes=None,
            report_genes=True
        )

    rank_hypergeometric_results = pd.concat(rank_hypergeometric_results, axis=0).reset_index().rename(
        {'level_0': 'Collection', 'level_1': 'Term'}, axis=1
    )
    rank_hypergeometric_results.rename({
        'Term': 'Geneset', 'n_overlap': 'OverlapSize', 'effective_set_size': 'EffectiveSetSize', 'original_set_size': 'OriginalSetSize', 'odds_ratio': 'OddsRatio', 'pval': 'FishersPValue', 'overlap': 'OverlappingGenes'
    }, axis=1).loc[
        :, [
            'Collection', 'Geneset', 'OverlapSize', 'EffectiveSetSize', 'OriginalSetSize', 'OddsRatio', 'FishersPValue', 'OverlappingGenes'
        ]
    ].to_csv(os.path.join(PROCESSED_DIR, 'organoids_vs_2d_gene_dependency_tail_enrichment_tests_full.csv'), index=False)
    
##### Minipool analysis #####
# setup for Chronos

def generate_sequence_map(sequence_metadata, pdna_sequence_ids=['HEK_CasNull_MP_control']):
    columns_to_tack = ['DepMap ID', 'Original Growth Pattern', 'Media condition', 'Serum status', 'Culture medium', 'Growth Format']
    if 'OncotreeLineage' in sequence_metadata.columns.tolist():
        columns_to_tack.append('OncotreeLineage')
    seq_map = sequence_metadata.loc[:, columns_to_tack]
    is_pdna = seq_map['Original Growth Pattern'] == "pDNA"
    non_pdnas = seq_map[~is_pdna].copy()
    pdnas = seq_map[is_pdna].loc[pdna_sequence_ids, :].copy()
    
    non_pdnas['sequence_ID'] = non_pdnas.index
    non_pdnas['cell_line_name'] = non_pdnas.index.map(lambda x: '_'.join(x.split('_')[:-1]))
    non_pdnas['pDNA_batch'] = 0
    non_pdnas['days'] = 21
    
    pdnas['sequence_ID'] = pdnas.index
    pdnas['cell_line_name'] = 'pDNA'
    pdnas['pDNA_batch'] = 0
    pdnas['days'] = 0
    
    primary_columns = ['sequence_ID', 'cell_line_name', 'pDNA_batch', 'days']
    secondary_columns = [c for c in seq_map.columns if c not in primary_columns]
    
    return pd.concat([non_pdnas, pdnas], axis=0).loc[:, primary_columns + secondary_columns]

## QC metrics

def compute_nnmd(df, poscons, negcons):
    pos = df.index.isin(poscons)
    neg = df.index.isin(negcons)
    nnmd = (df[pos].median() - df[neg].median()) / (df[neg] - df[neg].median()).abs().median()
    return nnmd

def compute_roc_auc(x, poscons, negcons):
    subset_srs = x[x.index.isin(poscons) | x.index.isin(negcons)].dropna()
    pos = (subset_srs.index.isin(negcons)).astype(int)
    if len(subset_srs) == 0 or len(pos) == 0:
        return np.nan
    fpr, tpr, _ = roc_curve(pos, subset_srs, pos_label=1)
    return auc(fpr, tpr)

# downstream analysis

def normalize_chronos_gene_effects(gene_effect_matrix, 
                                   negative_control_genes,
                                   positive_control_genes,
                                   individual=False):
    ge_matrix = gene_effect_matrix.copy()
    if individual:
        ge_matrix_T = ge_matrix.T
        ge_matrix_T -= ge_matrix_T.reindex(index=negative_control_genes).median()
        ge_matrix_T /= ge_matrix_T.reindex(index=positive_control_genes).abs().median()
        ge_matrix = ge_matrix_T.T
    else:
        ge_matrix -= ge_matrix.reindex(columns=negative_control_genes).median(axis=1).median()
        ge_matrix /= ge_matrix.reindex(columns=positive_control_genes).median(axis=1).abs().median()
    return ge_matrix

def gene_effect_linear_regression(gene_dependency_matrix, metadata_table, exogenous_variables, 
                                  combination_variables=[]):
    # format inputs so that patsy (in statsmodels) can recognize the formula
    original_variable_to_formatted_variable = {}
    for var in exogenous_variables:
        fixed_varname = var.replace(' ', '_')
        if metadata_table[var].dtype == np.dtype(float):
            original_variable_to_formatted_variable[var] = fixed_varname
        elif metadata_table[var].dtype == np.dtype('O'):
            original_variable_to_formatted_variable[var] = f'C({fixed_varname})'
    for combo in combination_variables:
        for var in combo:
            fixed_varname = var.replace(' ', '_')
            if metadata_table[var].dtype == np.dtype(float):
                original_variable_to_formatted_variable[var] = fixed_varname
            elif metadata_table[var].dtype == np.dtype('O'):
                original_variable_to_formatted_variable[var] = f'C({fixed_varname})'
    formatted_variable_to_original_variable = {v:k for k,v in original_variable_to_formatted_variable.items()}
    
    # single-variable terms
    independent_formula = ' + '.join([original_variable_to_formatted_variable[v] for v in exogenous_variables])
    
    # interaction_terms
    if len(combination_variables) > 0:
        interaction_terms = [' : '.join([original_variable_to_formatted_variable[v] for v in combo]) 
                             for combo in combination_variables]
        interaction_formula = ' + '.join(interaction_terms)
        base_formula = ' + '.join([independent_formula, interaction_formula])
    else:
        base_formula = independent_formula

    # fix the columns in the metadata_table to use the patsy-formatted names
    ## TODO: this has only been tested on categorical/str columns - don't want to strip the "C(",")" labels for numerics
    formatted_design_df = metadata_table.rename(
        {k:v[2:-1] for k,v in original_variable_to_formatted_variable.items()}, axis=1
    )
    
    # for each column in the gene_dependency_matrix, fit a linear model
    all_gene_results = dict()
    for gene in gene_dependency_matrix.columns:
        # fill the formula with the gene ("Q('')" causes patsy to recognize the input as a string, even with spaces/hyphens)
        formula = f"Q('{gene}') ~ {base_formula}"
        
        # extract the response variable to be joined with the predictors
        endog_var = gene_dependency_matrix.loc[:, gene]
        
        # fit the model
        linear_model = smf.glm(
            formula = formula,
            data = pd.concat([endog_var, formatted_design_df], axis=1)
        ).fit()
        
        # extract the coefficients and p-values from the model
        coefficients = linear_model.params.rename('coefficient')
        pvalues = linear_model.pvalues.rename('pvalue')
        all_gene_results[gene] = pd.concat([coefficients, pvalues], axis=1)
        
    all_gene_results = pd.concat(all_gene_results).reset_index().rename(
        {'level_0': 'gene', 'level_1': 'variable'}, axis=1
    )
    
    # return as a summary dict
    coefficients = all_gene_results.pivot(index='variable', columns='gene', values='coefficient')
    pvals = all_gene_results.pivot(index='variable', columns='gene', values='pvalue')
    return {'coefficients': coefficients, 'pvalues': pvals}

def f_test_variance(x, y, alternative='two-sided'):
    varx = np.var(x)
    vary = np.var(y)
    dfx = len(x) - 1
    dfy = len(y) - 1
    F = varx / vary
    pval = scipy.stats.f.cdf(F, dfx, dfy)
    if alternative == 'two-sided':
        pval = 2 * np.minimum(pval, 1-pval)
    elif alternative == 'greater':
        pval = 1 - pval
    elif alternative != 'less':
        raise ValueError('unrecognized test direction')
    return pval

# task functions

def run_minipool_qc(readcounts, sample_metadata, guide_metadata, test_guides, neg_ctrl_guides, pos_ctrl_guides):
    full_sequence_map = generate_sequence_map(sequence_metadata=sample_metadata, pdna_sequence_ids=['HEK_CasNull_MP_control']).reset_index(drop=True)
    full_guide_map = guide_metadata.rename({'Construct_Barcode': 'sgrna', 'Gene_ID': 'gene'}, axis=1)
    full_raw_readcounts = readcounts.loc[full_guide_map['sgrna'], full_sequence_map['sequence_ID']].T
    
    full_norm_readcounts = normalize_readcounts(
        full_raw_readcounts, 
        negative_control_sgrnas=neg_ctrl_guides, 
        sequence_map=full_sequence_map
    ).T

    full_sequence_guide_lfc = np.log2(
        calculate_fold_change(
            full_norm_readcounts.T,
            full_sequence_map,
            rpm_normalize=False
        )
    )

    full_screen_guide_lfc = full_sequence_guide_lfc.groupby(
        full_sequence_map.drop_duplicates(['sequence_ID']).set_index('sequence_ID')['cell_line_name']
    ).mean()

    replicate_map = full_sequence_map.groupby('cell_line_name').apply(lambda x: list(x['sequence_ID'])).drop('pDNA')
    replicate_dict = replicate_map[replicate_map.apply(lambda x: len(x) > 1)].to_dict()

    replicate_correlations_lfc_test_guides = dict()
    for screen, replicates in replicate_dict.items():
        replicate_correlations_lfc_test_guides[screen] = full_sequence_guide_lfc.loc[replicates, test_guides].T.corr().iloc[0, 1]
    replicate_correlations_lfc_test_guides = pd.Series(replicate_correlations_lfc_test_guides)
    
    sample_metadata['Total reads'] = readcounts.sum().reindex(index=sample_metadata.index)
    sample_metadata['NNMD'] = compute_nnmd(full_sequence_guide_lfc.T, poscons=pos_ctrl_guides, negcons=neg_ctrl_guides).reindex(sample_metadata.index)
    sample_metadata['ROC-AUC'] = full_sequence_guide_lfc.T.apply(lambda x: compute_roc_auc(x, poscons=pos_ctrl_guides, negcons=neg_ctrl_guides), axis=0).reindex(sample_metadata.index)
    sample_metadata['Replicate LFC Correlation'] = replicate_correlations_lfc_test_guides.reindex(
        index=full_sequence_map.set_index('sequence_ID').reindex(index=sample_metadata.index.tolist())['cell_line_name'].tolist()
    ).values.tolist()
    
    full_screen_metadata = full_sequence_map.drop_duplicates('cell_line_name').drop('sequence_ID', axis=1).set_index('cell_line_name').drop('pDNA')
    full_screen_metadata['Replicate LFC Correlation'] = replicate_correlations_lfc_test_guides.reindex(index=full_screen_metadata.index)

    full_screen_metadata['Total reads'] = readcounts.T.groupby(
        full_sequence_map.drop_duplicates(['sequence_ID']).set_index('sequence_ID')['cell_line_name']
    ).sum().sum(axis=1)

    full_screen_metadata['NNMD'] = compute_nnmd(full_screen_guide_lfc.T, poscons=pos_ctrl_guides, negcons=neg_ctrl_guides).reindex(full_screen_metadata.index)
    full_screen_metadata['ROC-AUC'] = full_screen_guide_lfc.T.apply(lambda x: compute_roc_auc(x, poscons=pos_ctrl_guides, negcons=neg_ctrl_guides), axis=0).reindex(full_screen_metadata.index)
    
    return sample_metadata, full_screen_metadata

def filter_samples(sample_metadata_with_qc):
    samples_to_keep =sample_metadata_with_qc[
        ((sample_metadata_with_qc['Total reads'] > 2.5e6) &
         (sample_metadata_with_qc['NNMD'] < -1.5) &
         ((sample_metadata_with_qc['Replicate LFC Correlation'] > 0.2) | 
          (sample_metadata_with_qc['Replicate LFC Correlation'].isna()))) |
        (sample_metadata_with_qc['Original Growth Pattern'] == "pDNA")
    ]

    replicate_groupings = samples_to_keep.groupby(['DepMap ID', 'Growth Format', 'Serum status']).count()

    cell_lines_with_serum_status_pairs = [x[0] for x in (replicate_groupings.groupby(
        level=[0, 1]).count().iloc[:, 0] > 1).loc[lambda x: x == True].index]
    cell_lines_with_growth_format_pairs = [x[0] for x in (replicate_groupings.groupby(
        level=[0, 2]).count().iloc[:, 0] > 1).loc[lambda x: x == True].index]

    samples_to_keep = samples_to_keep[
        samples_to_keep['DepMap ID'].isin(cell_lines_with_growth_format_pairs) |
        samples_to_keep['DepMap ID'].isin(cell_lines_with_serum_status_pairs) |
        (samples_to_keep['Original Growth Pattern'] == "pDNA")
    ]
    
    return samples_to_keep

def run_chronos(readcounts, sample_metadata_to_keep, guide_metadata, neg_ctrl_guides, neg_ctrl_genes, pos_ctrl_genes, seed=41):
    sequence_map = generate_sequence_map(
        pdna_sequence_ids=['HEK_CasNull_MP_control'],
        sequence_metadata=sample_metadata_to_keep
    )
    guide_map = guide_metadata.rename({'Construct_Barcode': 'sgrna', 'Gene': 'gene'}, axis=1)
    chronos_readcounts = readcounts.loc[guide_map['sgrna'], sequence_map['sequence_ID']].T
    screen_metadata = sequence_map.drop_duplicates('cell_line_name').drop('sequence_ID', axis=1).set_index('cell_line_name').drop('pDNA')
    
    np.random.seed(seed)
    model = Chronos(
        readcounts={'adhesome': readcounts.loc[:, sequence_map['sequence_ID']].T},
        guide_gene_map={'adhesome': guide_map},
        sequence_map={'adhesome': sequence_map},
        negative_control_sgrnas={'adhesome': neg_ctrl_guides},
        gene_effect_hierarchical=0.1, # default 0.1
        gene_effect_smoothing=3.0, # default 1.5
        kernel_width=5, # default 50
    )

    model.train(nepochs=1001)

    chronos_ge_unnormalized = model.gene_effect
    chronos_gene_effects = normalize_chronos_gene_effects(model.gene_effect, negative_control_genes=neg_ctrl_genes, positive_control_genes=pos_ctrl_genes)
    
    return chronos_gene_effects, screen_metadata, sequence_map
    
def run_linear_regression_analysis(chronos_gene_effects, full_screen_metadata, test_genes):
    linear_regression_results = gene_effect_linear_regression(
        chronos_gene_effects, 
        metadata_table=full_screen_metadata, 
        exogenous_variables=['Serum status', 'Growth Format', 'DepMap ID'],
        combination_variables=[('Growth Format', 'Serum status')],
    )

    growth_format_differences = pd.concat([
        linear_regression_results['coefficients'].loc["C(Growth_Format)[T.Plastic]"].rename('coefficient'),
        linear_regression_results['pvalues'].loc["C(Growth_Format)[T.Plastic]"].rename('pvalue')
    ], axis=1)
    growth_format_differences['coefficient'] = -growth_format_differences['coefficient']
    growth_format_differences = growth_format_differences.loc[test_genes]
    growth_format_differences['FDR'] = scipy.stats.false_discovery_control(growth_format_differences['pvalue'])

    serum_differences = pd.concat([
        linear_regression_results['coefficients'].loc["C(Serum_status)[T.Serum-free]"].rename('coefficient'),
        linear_regression_results['pvalues'].loc["C(Serum_status)[T.Serum-free]"].rename('pvalue')
    ], axis=1)
    serum_differences = serum_differences.loc[test_genes]
    serum_differences['FDR'] = scipy.stats.false_discovery_control(serum_differences['pvalue'])
    
    return linear_regression_results, growth_format_differences, serum_differences

def compute_geneset_variance_tests(growth_format_tests, serum_status_tests, test_gene_annotations):
    geneset_variance_df = {
        'Geneset':[], 
        'Size':[],
        'GrowthVariance': [],
        'GrowthOutVariance': [],
        'MediaVariance': [],
        'MediaOutVariance': [],
        'GrowthF': [],
        'MediaF': [],
        'GrowthPValue': [],
        'MediaPValue': [],
    }

    for gs in test_gene_annotations.columns.tolist():
        gs_genes = test_gene_annotations[test_gene_annotations[gs]].index.tolist()
        gs_out_genes = test_gene_annotations[~test_gene_annotations[gs]].index.tolist()
        growth_coefficients = growth_format_tests['coefficient'].reindex(index=gs_genes)
        growth_out_coefficients = growth_format_tests['coefficient'].reindex(index=gs_out_genes)
        serum_coefficients = serum_status_tests['coefficient'].reindex(index=gs_genes)
        serum_out_coefficients = serum_status_tests['coefficient'].reindex(index=gs_out_genes)

        geneset_variance_df['Geneset'].append(gs)
        geneset_variance_df['Size'].append(len(gs_genes))
        geneset_variance_df['GrowthVariance'].append(growth_coefficients.var())
        geneset_variance_df['GrowthOutVariance'].append(growth_out_coefficients.var())
        geneset_variance_df['MediaVariance'].append(serum_coefficients.var())
        geneset_variance_df['MediaOutVariance'].append(serum_out_coefficients.var())
        geneset_variance_df['GrowthF'].append(growth_coefficients.var() / growth_out_coefficients.var())
        geneset_variance_df['MediaF'].append(serum_coefficients.var() / serum_out_coefficients.var())
        geneset_variance_df['GrowthPValue'].append(f_test_variance(growth_coefficients, growth_out_coefficients, alternative='greater'))
        geneset_variance_df['MediaPValue'].append(f_test_variance(serum_coefficients, serum_out_coefficients, alternative='greater'))

    geneset_variance_df = pd.DataFrame.from_dict(geneset_variance_df)
    
    geneset_variance_df.loc[
        ~geneset_variance_df['GrowthPValue'].isna(),
        'GrowthFDR'
    ] = scipy.stats.false_discovery_control(geneset_variance_df['GrowthPValue'].dropna())

    geneset_variance_df.loc[
        ~geneset_variance_df['MediaPValue'].isna(),
        'MediaFDR'
    ] = scipy.stats.false_discovery_control(geneset_variance_df['MediaPValue'].dropna())
    
    return geneset_variance_df
    
def save_minipool_data(full_screen_metadata, mp_screen_metadata, chronos_gene_effect, linear_regression_results, growth_format_differences, serum_differences, geneset_variance_df, max_gene_effect_size=-0.5):
    mp_library_composition = load_file('minipool_library_composition.csv', index_col=0)

    chronos_gene_effect.to_csv(os.path.join(PROCESSED_DIR, 'minipool_chronos_gene_effects.csv'))
    
    any_dependency = (chronos_gene_effect < max_gene_effect_size).any().rename('AnyDependency')

    growth_format_differences.assign(
        DomeMean = chronos_gene_effect.loc[
            mp_screen_metadata.query('`Growth Format` == "Dome"').index
        ].mean().reindex(growth_format_differences.index),
        PlasticMean = chronos_gene_effect.loc[
            mp_screen_metadata.query('`Growth Format` == "Plastic"').index
        ].mean().reindex(growth_format_differences.index)
    ).join(
        any_dependency
    ).join(
        combine_annotations(
            mp_library_composition, annotation_columns=['Integrin', 'Actin Regulation', 'Adherens Junction', 'Tight Junction']
        ).rename('RepresentativeGeneset').fillna('Other')
    ).reset_index().rename(
        {'gene': 'Gene', 'coefficient': 'GrowthCoefficient', 'pvalue': 'GrowthPValue', 'FDR': 'GrowthFDR'}, axis=1
    ).to_csv(os.path.join(PROCESSED_DIR, 'minipool_growth_format_tests.csv'), index=False)

    serum_differences.assign(
        SerumMean = chronos_gene_effect.loc[
            mp_screen_metadata.query('`Serum status` == "Serum"').index
        ].mean().reindex(serum_differences.index),
        SerumFreeMean = chronos_gene_effect.loc[
            mp_screen_metadata.query('`Serum status` == "Serum-free"').index
        ].mean().reindex(serum_differences.index)
    ).join(
        any_dependency
    ).join(
        combine_annotations(
            mp_library_composition, annotation_columns=['Actin Regulation', 'Actin Nucleation', 'Adherens Junction', 'Tight Junction', 'Lipid Metabolism', 'Wnt Signaling']
        ).rename('RepresentativeGeneset').fillna('Other')
    ).reset_index().rename(
        {'gene': 'Gene', 'coefficient': 'MediaCoefficient', 'pvalue': 'MediaPValue', 'FDR': 'MediaFDR'}, axis=1
    ).to_csv(os.path.join(PROCESSED_DIR, 'minipool_serum_status_tests.csv'), index=False)
    
    geneset_variance_df.to_csv(os.path.join(PROCESSED_DIR, 'minipool_geneset_variance_tests.csv'), index=False)
    
def process_minipool_screens():
    readcounts = load_file('minipool_readcounts.csv', index_col=0)
    sample_metadata = load_file('minipool_sequence_metadata.csv', index_col=0)
    guide_metadata = load_file('minipool_guide_metadata.csv')
    test_gene_annotations = load_file('minipool_library_composition.csv', index_col=0).drop('Other', axis=1)
    
    pos_ctrl_guides = guide_metadata.query('Gene_category == "Positive_control"')['Construct_Barcode'].tolist()
    neg_ctrl_guides = guide_metadata.query('Gene_category == "Negative_control"')['Construct_Barcode'].tolist()
    pos_ctrl_genes = guide_metadata.query('Gene_category == "Positive_control"')['Gene'].unique().tolist()
    neg_ctrl_genes = guide_metadata.query('Gene_category == "Negative_control"')['Gene'].unique().tolist()
    test_guides = guide_metadata.query('Gene_category == "Test_gene"')['Construct_Barcode'].tolist()
    test_genes = guide_metadata.query('Gene_category == "Test_gene"')['Gene'].unique().tolist()
    
    sample_metadata_with_qc, full_screen_metadata = run_minipool_qc(readcounts, sample_metadata, guide_metadata, test_guides, neg_ctrl_guides, pos_ctrl_guides)
    
    sample_metadata_to_keep = filter_samples(sample_metadata_with_qc)
    
    chronos_gene_effects, mp_screen_metadata, sequence_map = run_chronos(readcounts, sample_metadata_to_keep, guide_metadata, neg_ctrl_guides, neg_ctrl_genes, pos_ctrl_genes, seed=CHRONOS_SEED)
    
    linear_regression_results, growth_format_differences, serum_differences = run_linear_regression_analysis(chronos_gene_effects, full_screen_metadata, test_genes)
    
    geneset_variance_df = compute_geneset_variance_tests(growth_format_differences, serum_differences, test_gene_annotations)
    
    save_minipool_data(full_screen_metadata, mp_screen_metadata, chronos_gene_effects, linear_regression_results, growth_format_differences, serum_differences, geneset_variance_df)

##### Genomic alterations #####

def configure_frequency_table(genes, lineages, alteration_types, tcga_sources, feature_labels):
    df_dict = {
        'Feature': [],
        'FeatureLabel': [],
        'FeatureID': [],
        'OncotreeLineage': [],
        'AlterationType': [],
        'TCGASource': []
    }
    for g, lins, alt, source, label in zip(genes, lineages, alteration_types, tcga_sources, feature_labels):
        if type(lins) != list:
            lins = [lins]
        for lin in lins:
            df_dict['Feature'].append(g)
            df_dict['OncotreeLineage'].append(lin)
            df_dict['AlterationType'].append(alt)
            df_dict['TCGASource'].append(source)
            df_dict['FeatureLabel'].append(label)
            
            if alt in ["HOT", "DAM"]:
                df_dict['FeatureID'].append(search_gene(g))
            else:
                df_dict['FeatureID'].append(g)
    return pd.DataFrame.from_dict(df_dict)

def compute_tcga_frequencies(alteration_metadata_df, 
                             sourced_data={
                                 ('TERT', 'CNS/Brain'): 0.8,
                                 ('18q', 'Bowel'): 0.66,
                                 ('9p', 'Pancreas'): 0.48,
                                 ('NRAS', 'Bowel'): 0.044,
                                 ('ESR1', 'Breast'): 0.07,
                                 ('SPOP', 'Prostate'): 0.13
                             }):
    supporting_data = load_file('tcga_select_mutations.csv')
    # for each (feature, lineage) pair
    any_mut_tcga_frequencies = dict()
    for i, r in alteration_metadata_df.iterrows():
        mut = r['Feature']
        lin = r['OncotreeLineage']
        if r['TCGASource'] == "Publication": # source the frequency from external data
            any_mut_tcga_frequencies[(mut, lin)] = sourced_data[(mut, lin)]
        else: # compute the frequency using a mutation matrix
            any_mut_tcga_frequencies[(mut, lin)] = (supporting_data[supporting_data['OncotreeLineage'] == lin][mut] != 'WT').mean()
            
    any_mut_tcga_frequencies = pd.Series(any_mut_tcga_frequencies).reset_index().rename({'level_0': 'Feature', 'level_1': 'OncotreeLineage', 0: 'Frequency'}, axis=1)
    return alteration_metadata_df.merge(any_mut_tcga_frequencies)

def compute_depmap_frequencies(alteration_metadata_df, model_subset, model_metadata, hotspot, damaging, arm_cn=None, 
                               supporting_data=None):
    # subset the corresponding mutation matrices
    alteration_type_to_feature_matrices = {
        "DAM": (damaging.reindex(index=model_subset, columns=alteration_metadata_df[alteration_metadata_df['AlterationType'] == "DAM"]['FeatureID'].unique().tolist()).dropna(how='all') > 0),
        "HOT": (hotspot.reindex(index=model_subset, columns=alteration_metadata_df[alteration_metadata_df['AlterationType'] == "HOT"]['FeatureID'].unique().tolist()).dropna(how='all') > 0),
    }
    
    # for each (feature, lineage) pair, use the appropriate matrix and compute the frequency within lineage
    depmap_frequencies = dict()
    for i, r in alteration_metadata_df.iterrows():
        mut = r['FeatureID']
        lin = r['OncotreeLineage']
        alt_type = r['AlterationType']
        if alt_type in ['Gain', 'Loss']:
            depmap_frequencies[(mut, lin)] = supporting_data[(mut, lin)]
        else:
            depmap_frequencies[(mut, lin)] = alteration_type_to_feature_matrices[r['AlterationType']].reindex(index=model_metadata[model_metadata['OncotreeLineage'] == lin].index.tolist())[mut].mean()
    depmap_frequencies = pd.Series(depmap_frequencies).reset_index().rename({'level_0': 'FeatureID', 'level_1': 'OncotreeLineage', 0: 'Frequency'}, axis=1)
    return alteration_metadata_df.merge(depmap_frequencies)

def compute_selected_alteration_frequencies():
    model_metadata = load_file('model_metadata.csv', index_col=0)
    hotspot = load_full_matrix('hotspot')
    damaging = load_full_matrix('damaging')
    
    frequency_metadata_df = configure_frequency_table(
        ['TP53', 'APC', 'TERT', 'KRAS', 'NRAS', 'ESR1', 'SPOP', '18q', '9p'],
        [['Bowel', 'Pancreas'], ['Bowel'], ['CNS/Brain'], ['Pancreas'], ['Bowel'], ['Breast'], ['Prostate'], ['Bowel'], ['Pancreas']],
        ['DAM', 'DAM', 'HOT', 'HOT', 'HOT', 'HOT', 'DAM', 'Loss', 'Loss'],
        ['cBioPortal', 'cBioPortal', 'Publication', 'cBioPortal', 'Publication', 'Publication', 'Publication', 'Publication', 'Publication'],
        ['TP53 damaging', 'APC damaging', 'TERT promoter hotspot', 'KRAS hotspot', 'NRAS hotspot', 'ESR1 hotspot', 'SPOP damaging', '18q loss', '9p loss']
    )

    alteration_heatmap_df = pd.concat([
        compute_tcga_frequencies(frequency_metadata_df).assign(Dataset='Tumor'),
        compute_depmap_frequencies(
            frequency_metadata_df, 
            model_metadata[model_metadata['IsNextGen']].index.tolist(),
            model_metadata, hotspot, damaging, arm_cn=None, supporting_data = {('18q', 'Bowel'): 21/25, ('9p', 'Pancreas'): 22/28}
        ).assign(TCGASource=np.nan, Dataset='DepMap NextGen'),
        compute_depmap_frequencies(
            frequency_metadata_df,
            model_metadata[model_metadata['IsTraditional2D']].index.tolist(),
            model_metadata, hotspot, damaging, arm_cn=None, supporting_data = {('18q', 'Bowel'): 58/69, ('9p', 'Pancreas'): 28/46}
        ).assign(TCGASource=np.nan, Dataset='DepMap Traditional')
    ], axis=0)
    
    alteration_heatmap_df.to_csv(os.path.join(PROCESSED_DIR, 'selected_genomic_alteration_frequencies.csv'), index=False)
        
##### Preprocess dependency tables #####

def preprocess_addictions(fdr_threshold=0.005, correlation_threshold=0.3):
    tda_summary = load_file('tda_dependency_classification.csv', local_dir=PROCESSED_DIR, index_col=0)
    geneset_table = load_file('geneset_table_updated.csv')
    next_gen_expr_addictions = load_file('expression_addiction_tests_full.csv', local_dir=PROCESSED_DIR, index_col=0)
    
    gene_functional_class_df = geneset_table[
        geneset_table['Geneset'].isin(['Oncogenes', 'Transcription factors'])
    ].assign(indicator=True).pivot(index='cds_gene_id', columns='Geneset', values='indicator').fillna(False)

    gene_functional_class_df.loc[(gene_functional_class_df['Oncogenes'] & gene_functional_class_df['Transcription factors']), 'FunctionalClass'] = 'Oncogenic TF'
    gene_functional_class_df.loc[(gene_functional_class_df['Oncogenes'] & ~gene_functional_class_df['Transcription factors']), 'FunctionalClass'] = 'Oncogene'
    gene_functional_class_df.loc[(~gene_functional_class_df['Oncogenes'] & gene_functional_class_df['Transcription factors']), 'FunctionalClass'] = 'TF'
    gene_functional_class_df.loc[(~gene_functional_class_df['Oncogenes'] & ~gene_functional_class_df['Transcription factors']), 'FunctionalClass'] = 'Other'
    
    next_gen_expr_addiction_volcano_df = next_gen_expr_addictions.join(tda_summary['NextGenClass']).dropna(subset=['NextGenPValue'])
    next_gen_expr_addiction_volcano_df = next_gen_expr_addiction_volcano_df[next_gen_expr_addiction_volcano_df['NextGenClass'].isin(['Strongly Selective', 'High Variance', 'Pan-dependency'])]
    next_gen_expr_addiction_volcano_df = next_gen_expr_addiction_volcano_df.join(gene_functional_class_df['FunctionalClass'])
    next_gen_expr_addiction_volcano_df['FunctionalClass'] = next_gen_expr_addiction_volcano_df['FunctionalClass'].fillna('Other')
    next_gen_expr_addiction_volcano_df['FDR'] = scipy.stats.false_discovery_control(next_gen_expr_addiction_volcano_df['NextGenPValue'])
    next_gen_expr_addiction_volcano_df['Significant'] = (
        (next_gen_expr_addiction_volcano_df['FDR'] < fdr_threshold) &
        (next_gen_expr_addiction_volcano_df['NextGenPearsonR'] < -correlation_threshold)
    )
    
    next_gen_expr_addiction_volcano_df.to_csv(os.path.join(PROCESSED_DIR, 'expression_addiction_table.csv'))
    
    return next_gen_expr_addiction_volcano_df

def preprocess_paralog_associations(variance_threshold=1, fdr_threshold=0.005, correlation_threshold=0.3):
    tda_summary = load_file('tda_dependency_classification.csv', local_dir=PROCESSED_DIR, index_col=0)
    next_gen_paralogs = load_file('paralog_dependency_tests_full.csv', local_dir=PROCESSED_DIR)
    next_gen_paralog_dependency_volcano_df = next_gen_paralogs.merge(tda_summary['NextGenClass'], left_on='Dependency', right_index=True, how='left')
    next_gen_paralog_dependency_volcano_df = next_gen_paralog_dependency_volcano_df[next_gen_paralog_dependency_volcano_df['NextGenClass'].isin(['Strongly Selective', 'High Variance', 'Pan-dependency'])]
    next_gen_paralog_dependency_volcano_df = next_gen_paralog_dependency_volcano_df[next_gen_paralog_dependency_volcano_df['BiomarkerVariance'] > variance_threshold]
    next_gen_paralog_dependency_volcano_df['FDR'] = scipy.stats.false_discovery_control(next_gen_paralog_dependency_volcano_df['NextGenPValue'])
    next_gen_paralog_dependency_volcano_df['Significant'] = (
        (next_gen_paralog_dependency_volcano_df['NextGenPearsonR'] > correlation_threshold) &
        (next_gen_paralog_dependency_volcano_df['FDR'] < fdr_threshold)
    )
    
    next_gen_paralog_dependency_volcano_df.to_csv(os.path.join(PROCESSED_DIR, 'paralog_dependency_table.csv'), index=False)
    
    return next_gen_paralog_dependency_volcano_df

def summarize_lineage_test_range(lineage_test_df, groupby_column='Dependency', stat_column='Mean Difference', additional_columns=['OncotreeLineage', 'Ingroup n'], require_min_ingroup_n=0):
    max_idxs = lineage_test_df[lineage_test_df['Ingroup n'] >= require_min_ingroup_n].groupby(groupby_column)[stat_column].idxmax()
    max_stats = lineage_test_df.loc[max_idxs.values.tolist(), [groupby_column, stat_column] + additional_columns].rename(
        {c: 'Max ' + c for c in [stat_column] + additional_columns}, axis=1
    )
    
    min_idxs = lineage_test_df[lineage_test_df['Ingroup n'] >= require_min_ingroup_n].groupby(groupby_column)[stat_column].idxmin()
    min_stats = lineage_test_df.loc[min_idxs.values.tolist(), [groupby_column, stat_column] + additional_columns].rename(
        {c: 'Min ' + c for c in [stat_column] + additional_columns}, axis=1
    )
    
    summary_stat_df = max_stats.merge(min_stats, left_on=groupby_column, right_on=groupby_column)
    return summary_stat_df

def preprocess_onc_tsg_associations(fdr_threshold=0.05, mean_difference_threshold=0.3, min_lineage_size=5):
    cancer_gene_annotations = load_file('next_gen_cancer_gene_annotations.csv', local_dir=PROCESSED_DIR, index_col=0)
    tda_summary = load_file('tda_dependency_classification.csv', local_dir=PROCESSED_DIR, index_col=0)
    next_gen_onc_tsg_tests = load_file('onc_tsg_alteration_dependency_tests_full.csv', local_dir=PROCESSED_DIR)
    next_gen_dependency_lineage_tests = load_file('next_gen_dependency_lineage_tests.csv', local_dir=PROCESSED_DIR)

    next_gen_onc_tsg_dependency_volcano_df = next_gen_onc_tsg_tests.merge(
        tda_summary['NextGenClass'], left_on='Dependency', right_index=True, how='left'
    ).merge(
        cancer_gene_annotations['AnyHotspotInNextGen'], left_on='Biomarker', right_index=True
    )
    next_gen_onc_tsg_dependency_volcano_df = next_gen_onc_tsg_dependency_volcano_df[
        next_gen_onc_tsg_dependency_volcano_df['NextGenClass'].isin(['Strongly Selective', 'High Variance', 'Pan-dependency']) & 
        next_gen_onc_tsg_dependency_volcano_df['AnyHotspotInNextGen']
    ]
    next_gen_onc_tsg_dependency_volcano_df = next_gen_onc_tsg_dependency_volcano_df.merge(
        summarize_lineage_test_range(next_gen_dependency_lineage_tests, require_min_ingroup_n=min_lineage_size).rename({
            'Max Mean Difference':'MaxLineageMeanDifference', 'Max OncotreeLineage':'MaxLineage', 'Min Mean Difference':'MinLineageMeanDifference', 'Min OncotreeLineage':'MinLineage'
        }, axis=1).loc[:, ['Dependency', 'MaxLineageMeanDifference', 'MaxLineage', 'MinLineageMeanDifference', 'MinLineage']]
    )
    next_gen_onc_tsg_dependency_volcano_df['FDR'] = scipy.stats.false_discovery_control(next_gen_onc_tsg_dependency_volcano_df['NextGenMannWhitneyPValue'])
    next_gen_onc_tsg_dependency_volcano_df['HitClass'] = 'Insignificant'
    next_gen_onc_tsg_dependency_volcano_df.loc[
        ((next_gen_onc_tsg_dependency_volcano_df['FDR'] < fdr_threshold) & (next_gen_onc_tsg_dependency_volcano_df['NextGenMeanDifference'] < -mean_difference_threshold)) &
        (next_gen_onc_tsg_dependency_volcano_df['NextGenMeanDifference'] < next_gen_onc_tsg_dependency_volcano_df['MinLineageMeanDifference']), 
        'HitClass'
    ] = 'Driven by alteration'
    next_gen_onc_tsg_dependency_volcano_df.loc[
        ((next_gen_onc_tsg_dependency_volcano_df['FDR'] < fdr_threshold) & (next_gen_onc_tsg_dependency_volcano_df['NextGenMeanDifference'] < -mean_difference_threshold)) &
        (next_gen_onc_tsg_dependency_volcano_df['NextGenMeanDifference'] > next_gen_onc_tsg_dependency_volcano_df['MinLineageMeanDifference']),
        'HitClass'
    ] = 'Driven by lineage'
    next_gen_onc_tsg_dependency_volcano_df.to_csv(os.path.join(PROCESSED_DIR, 'onc_tsg_dependency_table.csv'), index=False)
    
def preprocess_pdac_dependency_associations():
    tda_summary = load_file('tda_dependency_classification.csv', local_dir=PROCESSED_DIR, index_col=0)
    geneset_table = load_file('geneset_table_updated.csv')
    wnt_pathway = geneset_table[geneset_table['Geneset'] == 'GOBP_WNT_SIGNALING_PATHWAY']['cds_gene_id'].tolist()
    pdac_classical_associations = load_file('pdac_classical_dependency_tests_full.csv', local_dir=PROCESSED_DIR, index_col=0)
    pdac_classical_volcano_df = pdac_classical_associations.join(tda_summary['OrganoidsClass'])
    pdac_classical_volcano_df = pdac_classical_volcano_df[pdac_classical_volcano_df['OrganoidsClass'].isin(['Strongly Selective', 'High Variance', 'Pan-dependency'])]
    pdac_classical_volcano_df['FDR'] = scipy.stats.false_discovery_control(pdac_classical_volcano_df['NextGenPValue'])
    pdac_classical_volcano_df['Wnt Signaling'] = pdac_classical_volcano_df.index.isin(wnt_pathway)
    
    enrichment_x = pdac_classical_volcano_df['NextGenPearsonR'].sort_values()
    gsea_res = gseapy.prerank(rnk=enrichment_x, gene_sets={'Wnt Signaling': wnt_pathway})
    enrichment_y = [-x for x in gsea_res.results['Wnt Signaling']['RES'][::-1]]
    enrichment_gene_order = enrichment_x.index.tolist()
    
    pdac_classical_volcano_df = pdac_classical_volcano_df.join(pd.Series(enrichment_y, index=enrichment_gene_order, name='RunningEnrichmentScore'), )
    
    pdac_classical_volcano_df.to_csv(os.path.join(PROCESSED_DIR, 'pdac_classical_dependency_table.csv'))
    
def prepare_wnt_example_df(wnt_genes_to_regress):
    screen_gene_effect = load_file('screen_gene_effect.csv', index_col=0)
    screen_metadata = load_file('screen_metadata.csv', index_col=0)
    geneset_scores = load_file('celligner_geneset_scores.csv', local_dir=PROCESSED_DIR, index_col=0)
    pdac_classical_scores = geneset_scores['PDAC-classical']
    screen_metadata = screen_metadata.merge(pdac_classical_scores, left_on='ModelID', right_index=True).replace(lineage_replacement)

    organoid_screen_metadata = screen_metadata[(screen_metadata['IsNextGen']) & (screen_metadata['ScreenType'].isin(["2DO", "3DO"]))]
    organoid_lineages = organoid_screen_metadata['OncotreeLineage'].value_counts().index.tolist()
    adherent_organoid_lineage_matched_screen_metadata = screen_metadata[
        (screen_metadata['IsAdherent2D']) & 
        (screen_metadata['OncotreeLineage']).isin(organoid_lineages)
    ]

    # format the dataframe for plotting
    wnt_dependency_df = []
    for g in wnt_genes_to_regress:
        wnt_ex_df = pd.concat([
            pd.Series(g, index=screen_gene_effect.index.tolist(), name='Dependency'),
            screen_gene_effect.loc[:, g].rename('DependencyGeneEffect'),
            screen_metadata[['PDAC-classical', 'OncotreeLineage']],
            pd.concat([
                pd.Series('DepMap Organoid', index=organoid_screen_metadata.index.tolist(), name='SampleSet'),
                pd.Series('DepMap Adherent', index=adherent_organoid_lineage_matched_screen_metadata.index.tolist(), name='SampleSet')
            ], axis=0)
        ], axis=1)
        wnt_dependency_df.append(wnt_ex_df)
    wnt_dependency_df = pd.concat(wnt_dependency_df, axis=0).dropna(subset=['SampleSet'])
    wnt_dependency_df.rename_axis('ScreenID').to_csv(os.path.join(PROCESSED_DIR, 'wnt_example_table.csv'))
    
##### Global growth format analysis #####
    
def filter_geneset_results(ranked_diff_hypergeometric_results, min_size=50, max_size=250, fdr_threshold=0.05, top_n=None):
    significant_ranked_hypergeometric_results_df = ranked_diff_hypergeometric_results[
        (ranked_diff_hypergeometric_results['OriginalSetSize'] >= min_size) & 
        (ranked_diff_hypergeometric_results['OriginalSetSize'] <= max_size) &
        (ranked_diff_hypergeometric_results['EffectiveSetSize'] > 100) 
        # &(ranked_diff_hypergeometric_results['EffectiveSetSize'] <= max_size)
    ].copy()
    # recalculate FDRs after filtering
    significant_ranked_hypergeometric_results_df['FDR'] = scipy.stats.false_discovery_control(
        significant_ranked_hypergeometric_results_df['FishersPValue']
    )
    significant_ranked_hypergeometric_results_df.to_csv("tmp-files/geneset_w_fdr_200genes.csv")
    significant_ranked_hypergeometric_results_df = significant_ranked_hypergeometric_results_df[
        significant_ranked_hypergeometric_results_df['FDR'] < fdr_threshold
    ]
    significant_ranked_hypergeometric_results_df = significant_ranked_hypergeometric_results_df.sort_values(['FDR', 'OverlapSize'], ascending=[True, False])
    
    return significant_ranked_hypergeometric_results_df.head(top_n)

def cluster_genesets(significant_ranked_hypergeometric_results_df, combined_geneset_dict, min_similarity=0.3):
    # compute similarity using the overlap coefficient
    sig_geneset_similarity = geneset_min_coverage_similarity(combined_geneset_dict, geneset_subset=significant_ranked_hypergeometric_results_df['Geneset'].tolist())
    geneset_clusters = partition_by_connected_component(sig_geneset_similarity > min_similarity)
    # plot a heatmap for visualization
    cluster_sizes = geneset_clusters.value_counts()
    sorter_intra_mtx = sig_geneset_similarity.copy()
    sorter_inter_mtx = sig_geneset_similarity.copy()
    for i in range(sorter_intra_mtx.shape[0]):
        sorter_intra_mtx.iloc[i, i] = np.nan
        sorter_inter_mtx.iloc[i, i] = np.nan
        for j in range(sorter_intra_mtx.shape[0]):
            if sorter_intra_mtx.iloc[i, j] < min_similarity:
                sorter_intra_mtx.iloc[i, j] = np.nan
            else:
                sorter_inter_mtx.iloc[i, j] = np.nan
    intracluster_means = sorter_intra_mtx.mean()
    intercluster_means = sorter_inter_mtx.mean()
    sorter = pd.concat([
        geneset_clusters.replace(cluster_sizes).rename('cluster_size'), 
        intracluster_means.rename('intracluster_mean_similarity'),
        intercluster_means.rename('intercluster_mean_similarity'),
    ], axis=1)
    order = sorter.sort_values(['cluster_size', 'intracluster_mean_similarity', 'intercluster_mean_similarity'], ascending=[False, False, False]).index.tolist()
    sig_geneset_similarity.loc[order, order].to_csv(os.path.join(PROCESSED_DIR, 'organoids_vs_2d_top_geneset_similarity_matrix.csv'))
    
    return geneset_clusters
    
def preprocess_global_format_comparison():
    org_2d_diff_df = load_file('organoids_vs_2d_gene_dependency_tests_full.csv', local_dir=PROCESSED_DIR, index_col=0)
    org_2d_diff_df['GeneSymbol'] = org_2d_diff_df.index.str.split(' ').str[0]
    org_2d_diff_df = org_2d_diff_df[~org_2d_diff_df['DropEssential']]
    org_2d_diff_df['MannWhitneyFDR'] = scipy.stats.false_discovery_control(org_2d_diff_df['MannWhitneyPValue'])
    ranked_diff_hypergeometric_results = load_file('organoids_vs_2d_gene_dependency_tail_enrichment_tests_full.csv', local_dir=PROCESSED_DIR)
    significant_ranked_hypergeometric_results_df = filter_geneset_results(ranked_diff_hypergeometric_results)

    geneset_table = load_file('geneset_table_updated.csv')
    combined_geneset_dict = geneset_table[geneset_table['Source'].isin(['Hallmark', 'KEGG'])].groupby('Geneset').apply(lambda x: list(x['GeneSymbol'])).to_dict()
    
    geneset_clusters = cluster_genesets(significant_ranked_hypergeometric_results_df, combined_geneset_dict, min_similarity=0.3)
    org_2d_diff_geneset_df = significant_ranked_hypergeometric_results_df.merge(
        geneset_clusters.rename('Cluster'), left_on='Geneset', right_index=True
    )

    # create geneset priority order by significance
    target_geneset_order = []
    for cl in org_2d_diff_geneset_df['Cluster'].drop_duplicates().values.tolist():
        target_geneset_order.append(geneset_clusters.loc[geneset_clusters == cl].index.tolist()[0])
    
    # assign geneset annotations by cluster
    genes_to_clusters_df = pd.DataFrame(index=org_2d_diff_df['GeneSymbol'])
    for cl in org_2d_diff_geneset_df['Cluster'].unique().tolist():
        all_genes_in_cluster = set()
        genesets = org_2d_diff_geneset_df[org_2d_diff_geneset_df['Cluster'] == cl]['Geneset'].tolist()
        for gs in genesets:
            all_genes_in_cluster = all_genes_in_cluster.union(combined_geneset_dict[gs])
        genes_to_clusters_df[f'{cl}'] = genes_to_clusters_df.index.isin(all_genes_in_cluster)
    
    # prepare dataframe for plotting
    org_2d_diff_df = org_2d_diff_df.merge(
        combine_annotations(
            genes_to_clusters_df, annotation_columns=[str(geneset_clusters.loc[gs]) for gs in target_geneset_order]
        ).fillna('Other').rename('Cluster'),
        left_on = 'GeneSymbol', right_index=True
    )
    
    org_2d_diff_geneset_df.to_csv(os.path.join(PROCESSED_DIR, 'organoids_vs_2d_geneset_enrichment_table.csv'), index=False)
    org_2d_diff_df.to_csv(os.path.join(PROCESSED_DIR, 'organoids_vs_2d_dependency_volcano_table.csv'))
    
    
def main():
    # all major analysis: input -> results
    print('calling dependency classes...')
    call_gene_classes()
    print('computing geneset scores...')
    compute_geneset_scores()
    print('computing metaprogram tests...')
    compare_metaprograms()
    print('calculating celligner results...')
    run_celligner_classification()
    print('computing expression addiction tests...')
    run_expression_addiction_analysis()
    print('computing paralog dependency tests...')
    run_paralog_analysis()
    print('computing oncogene / tumor suppressor dependency tests...')
    run_oncogene_tsg_association_analysis()
    print('computing CDK6 biomarker tests...')
    compute_gbm_cdk6_dependency_biomarker_tests(correlation_threshold=0.3, correlation_significance_threshold=0.005, min_mutations=5, mean_effect_threshold=0.3, mwu_significance_threshold=0.05)
    print('computing PDAC-classical dependency tests...')
    run_pdac_classical_analysis()
    print('comparing organoid and adherent dependency scores...')
    compare_next_gen_vs_adherent_dependency()
    print('processing minipool screens...')
    process_minipool_screens()
    
    # clean up some files for convenience with figure generation
    # figure 1
    print('Figure 1: computing alteration frequencies across datasets...')
    compute_selected_alteration_frequencies()
    
    # figure 2
    print('Figure 2/4: pruning dependency tests...')
    preprocess_addictions(fdr_threshold=0.005, correlation_threshold=0.3)
    preprocess_paralog_associations(variance_threshold=1, fdr_threshold=0.005, correlation_threshold=0.3)
    preprocess_onc_tsg_associations(fdr_threshold=0.05, mean_difference_threshold=0.3, min_lineage_size=5)
    
    # figure 4
    preprocess_pdac_dependency_associations()
    prepare_wnt_example_df(['TCF7L2 (6934)', 'WLS (79971)', 'FZD5 (7855)',  'MESD (23184)'])
    
    # figure 5
    print('Figure 5: summarizing organoid vs adherent analysis...')
    preprocess_global_format_comparison()

if __name__ == "__main__":
    main()
    