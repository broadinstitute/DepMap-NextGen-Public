import os
import pandas as pd
import numpy as np
import requests
import urllib.request

import scipy.stats
import statsmodels.api as sm
import statsmodels.formula.api as smf

from tqdm import tqdm

from constants import *

##### File ingestion #####

FIGSHARE_AUTH_TOKEN = pd.read_csv(os.path.join(ASSETS_DIR, 'figshare_token.txt')).columns.tolist()[0] # for programmatic access to private figshare data
api_call_headers = {'Authorization': f'token {FIGSHARE_AUTH_TOKEN}'}

FIGSHARE_FILE_MAP = {
    ### Should be loaded into the PROCESSED_DIR folder ###
    #"filename": "token",
  
}

FULL_MATRIX_MAP = {
    'expression': ['traditional_expression.csv', 'next_gen_expression.csv'],
    'copy_number': ['traditional_copy_number.csv', 'next_gen_copy_number.csv'],
    'damaging': ['traditional_damaging.csv', 'next_gen_damaging.csv'],
    'hotspot': ['traditional_hotspot.csv', 'next_gen_hotspot.csv'],
    'omics_signatures': ['traditional_omics_signatures.csv', 'next_gen_omics_signatures.csv']
}

def download_file_from_figshare(filename, local_dir=INPUT_DIR):
    '''
    Download a file from figshare into the local data folder
    '''
    html = "https://ndownloader.figshare.com/files/"
    file_num = FIGSHARE_FILE_MAP[filename]
    r = requests.get(f'https://ndownloader.figshare.com/files/{file_num}', allow_redirects=True, headers=api_call_headers)
    
    if r.status_code != 200:
        print('Error:', r.content)
    else:
        print(f"Downloading {filename} from {html}{file_num}...")
        with open(os.path.join(local_dir, filename), 'wb') as f:
            for chunk in r.iter_content(1024):
                f.write(chunk)
        print(f'Saved to {os.path.join(local_dir, filename)}')
        
    return os.path.join(local_dir, filename)

def load_file(filename, local_dir=INPUT_DIR, **kwargs):
    '''
    Read a file in from a local data folder, and pull it from figshare if absent
    '''
    if filename not in os.listdir(local_dir):
        filepath = download_file_from_figshare(filename, local_dir=local_dir)
    else:
        filepath = os.path.join(local_dir, filename)
    table = pd.read_csv(filepath, **kwargs)
    return table

def load_full_matrix(full_matrix_name):
    '''
    For data that crosses DepMap datasets, read both files, align the features, and join the matrices
    '''
    model_metadata = load_file('model_metadata.csv').set_index('ModelID')
    matrices_to_stitch = []
    for f in FULL_MATRIX_MAP[full_matrix_name]:
        partial_matrix = load_file(f, index_col=0)
        if full_matrix_name in ['expression', 'copy_number', 'damaging', 'hotspot']:
            partial_matrix_rename = remap_genes(partial_matrix.columns.tolist(), allow_duplicated_out=False)
            partial_matrix = partial_matrix.reindex(columns=partial_matrix_rename.index.tolist()).rename(partial_matrix_rename, axis=1)
        matrices_to_stitch.append(partial_matrix)
    full_matrix = pd.concat(matrices_to_stitch, axis=0)
    full_matrix = full_matrix[~full_matrix.index.duplicated(keep='last')]
    full_matrix = full_matrix.reindex(index=model_metadata.index.tolist()).dropna(how='all', axis=0)
    return full_matrix

##### Gene manipulation #####

hgnc_annotations = load_file('hgnc_table.csv')
hgnc_annotations = hgnc_annotations[~hgnc_annotations["cds_gene_id"].isna()]
hgnc_annotations["cds_gene_id"] = hgnc_annotations["cds_gene_id"].astype(str)

# function to resolve gene identity crises
def remap_genes(gene_tag_list, global_id_mapping=hgnc_annotations[['symbol', 'entrez_id', 'alias_symbol', 'prev_symbol', 'cds_gene_id']], 
                id_column='entrez_id', symbol_column='symbol', output_column='cds_gene_id', alternative_symbol_columns=['alias_symbol', 'prev_symbol'],
                allow_duplicated_out=False, preserve_unmapped=True):
    '''
    Ingests a list of CDS-formatted gene tags, expected to be "<HGNC Symbol> (<Entrez ID>)", but will also 
    accept only symbols. Then maps the genes to the current stable names according to the provided gene table. 
    Matches genes in order of Entrez ID, followed by HGNC symbol, then by any of the alternative symbols.
    '''
    
    assert output_column in global_id_mapping.columns.tolist(), f'"`{output_column}`" not found in `global_id_mapping`, please specify an output format'
    global_id_mapping = global_id_mapping.copy()
    
    # if the desired output is one of the inputs, this resolves the problem of duplicated columns
    if output_column in [id_column, symbol_column]:
        _output_column = '_' + output_column
        global_id_mapping[_output_column] = global_id_mapping[output_column]
    else:
        _output_column = output_column
    
    # process the ground truth labels and summarize
    ground_truth_mapping = []
    print('Processing the ground truth mapping...')
    print(f'Aligning `{output_column}` using `{id_column}` and `{symbol_column}`...')
    ground_truth_mapping.append(global_id_mapping[[symbol_column, id_column, _output_column]].dropna())
    for alt_sym_col in alternative_symbol_columns:
        if alt_sym_col in global_id_mapping.columns.tolist():
            print(f'Supplementing additional symbols from `{alt_sym_col}`...')
            map_supplement = global_id_mapping[[alt_sym_col, id_column, _output_column]].dropna().rename({alt_sym_col: symbol_column}, axis=1)
            map_supplement[symbol_column] = map_supplement[symbol_column].str.split('|')
            map_supplement = map_supplement.explode(symbol_column)
            ground_truth_mapping.append(map_supplement)
    ground_truth_mapping = pd.concat(ground_truth_mapping, axis=0, ignore_index=True).drop_duplicates(subset=[symbol_column], keep='first')
    print(f'Found {len(ground_truth_mapping[symbol_column].unique())} unique entries for `{symbol_column}` mapping to {len(ground_truth_mapping[id_column].unique())} unique entries for `{id_column}`')
    
    # for the inputs, split on whitespace and interpret the components as symbol and id
    print('Processing the inputs...')
    mapping_to_fix = pd.DataFrame({
        'input_tag': gene_tag_list,
        'input_symbol': [x.split(' ')[0] for x in gene_tag_list],
        'input_id': [x.split(' ')[1].strip('()') if ' ' in x else 'Unknown' for x in gene_tag_list]
    })
    mapping_to_fix['input_id'] = mapping_to_fix['input_id'].replace({'Unknown': np.nan}).astype(float)
    
    output_mappings = []

    to_fix_by_id = mapping_to_fix[~mapping_to_fix['input_id'].isna()]
    to_fix_by_symbol = []
    to_fix_by_symbol.append(mapping_to_fix[mapping_to_fix['input_id'].isna()])
    
    # attempt to resolve genes with ids first
    if len(to_fix_by_id) > 0:
        to_fix_by_id = to_fix_by_id.merge(ground_truth_mapping[[id_column, _output_column]].drop_duplicates(), 
                                          left_on='input_id', right_on=id_column, how='left').drop(id_column, axis=1)
        mapped_by_id = to_fix_by_id[~to_fix_by_id[_output_column].isna()]
        output_mappings.append(mapped_by_id)
        print(f'Mapped {(len(mapped_by_id) / len(mapping_to_fix)) * 100 :.02f}% of inputs by `{id_column}`...')
        
        # send the genes that failed id matching to the symbol matching inputs
        unmapped = to_fix_by_id[to_fix_by_id[_output_column].isna()].drop(_output_column, axis=1)
        if len(unmapped) > 0:
            to_fix_by_symbol.append(unmapped)
            print(f'Attempting to map the remainder by `{symbol_column}`...')

    # resolve the remaining genes with symbols
    if len(to_fix_by_symbol) > 0:
        to_fix_by_symbol = pd.concat(to_fix_by_symbol, axis=0)
        to_fix_by_symbol = to_fix_by_symbol.merge(ground_truth_mapping[[symbol_column, _output_column]].drop_duplicates(), 
                                                  left_on='input_symbol', right_on=symbol_column, how='left').drop(symbol_column, axis=1)
        mapped_by_symbol = to_fix_by_symbol[~to_fix_by_symbol[_output_column].isna()]
        output_mappings.append(mapped_by_symbol)
        print(f'Mapped {(len(mapped_by_symbol) / len(mapping_to_fix)) * 100 :.02f}% of inputs by `{symbol_column}`...')
        
        # remaining genes have failed all mapping attempts, summarize below
        unmapped = to_fix_by_symbol[to_fix_by_symbol[_output_column].isna()]
        
    print(f'{(len(unmapped) / len(mapping_to_fix)) * 100 :.02f}% of inputs remain unmapped')
    if len(unmapped) > 0:
        print('\t Unmapped inputs: {}'.format(', '.join(unmapped['input_tag'].tolist())))
    
    output_mappings = pd.concat(output_mappings, axis=0, ignore_index=True)
    if not allow_duplicated_out: # keep only the first instance of duplicate output ids
        output_mappings = output_mappings.drop_duplicates(subset=_output_column, keep='first')
    if preserve_unmapped: # for unmapped genes, keep them in the output mapping in their input forms
        output_mappings = pd.concat([
            output_mappings,
            pd.DataFrame({'input_tag': unmapped['input_tag'].tolist(), _output_column: unmapped['input_tag'].tolist()})
        ], axis=0)
    
    return output_mappings.set_index('input_tag')[_output_column]

# helper function for quick searching genes
def search_gene(query_gene, gene_list=hgnc_annotations['cds_gene_id'].tolist(), prefix=None):
    if prefix is not None:
        search_result = [x for x in gene_list if query_gene in x and x.startswith(prefix)]
        exact_filter = [x for x in search_result if x.split(' ')[0][-len(query_gene):] == query_gene]
        if len(exact_filter) == 1:
            return exact_filter[0]
    else:
        search_result = [x for x in gene_list if query_gene in x]
        exact_filter = [x for x in search_result if x.split(' ')[0] == query_gene]
        if len(exact_filter) == 1:
            return exact_filter[0]
    if len(search_result) < 1:
        raise ValueError(f'No gene found: {query_gene}')
    if len(search_result) > 1:
        print('More than one gene found:', str(search_result), f'... using {search_result[0]}')
    return search_result[0]
    
# helper function for annotating loci
def add_cytobands(df, column_names=[], locus_mapping=hgnc_annotations.set_index('cds_gene_id')['location_sortable']):
    df = df.copy()
    for c in column_names:
        if c in df.columns.tolist():
            df[f'{c} Cytoband'] = df[c].map(locus_mapping.get)
    return df

def expand_model_matrix_to_screens(matrix, model_mapping):
    return matrix.merge(
        model_mapping.to_frame(), how='inner', left_index=True, right_index=True
    ).set_index(model_mapping.name, drop=True)


##### Statistical testing #####

# correlation functions
def group_cols_with_same_mask(x):
    """
    Group columns with the same indexes of NAN values.
    
    Return a sequence of tuples (mask, columns) where columns are the column indices
    in x which all have the mask.
    """
    per_mask = {}
    for i in range(x.shape[1]):
        o_mask = np.isfinite(x[:, i])
        o_mask_b = np.packbits(o_mask).tobytes()
        if o_mask_b not in per_mask:
            per_mask[o_mask_b] = [o_mask, []]
        per_mask[o_mask_b][1].append(i)
    return per_mask.values()

def quick_corr(x, y, n=10):
    assert isinstance(x, pd.Series), "x must be a series"
    assert isinstance(y, pd.DataFrame), "y must be a df"
    x, y = x.align(y, join="outer")
    x_mask = ~x.isna().values
    y_groups = group_cols_with_same_mask(y.values)
    all_corrs = []
    for y_mask, y_columns in y_groups:
        combined_mask = x_mask & y_mask
        x_aligned = x[combined_mask]
        y_aligned = y.loc[combined_mask, :].iloc[:, y_columns]
        
        if y.shape[1] == 1:
            crs = scipy.stats.pearsonr(x_aligned.values, y_aligned.iloc[:, 0])
            corrs = pd.DataFrame(crs, index=['Pearson Correlation', 'p-value'], columns=y_aligned.columns.tolist()).T
        else:
            corrs = y_aligned.apply(lambda z: scipy.stats.pearsonr(x_aligned, z)).T.rename(columns={0:"Pearson Correlation", 1:'p-value'}).dropna()
        
        corrs['n'] = np.sum(combined_mask)
        all_corrs.append(corrs)
    all_corrs = pd.concat(all_corrs, axis=0)
    return all_corrs.sort_values('p-value').head(n)

def get_top_correlation_biomarkers(gene, 
                                   sample_subset,
                                   gene_effect_matrix, 
                                   biomarker_matrix,
                                   ntop=10):
    precorr = fast_cor(
        gene_effect_matrix.loc[sample_subset, [gene]], 
        biomarker_matrix
    ).T
    top_corrs = precorr.loc[:, gene].abs().sort_values(ascending=False).head(ntop).index.tolist()
    
    corrs = quick_corr(
        gene_effect_matrix.loc[sample_subset, gene].dropna(), biomarker_matrix.loc[:, top_corrs], 
        n=ntop
    )
    
    corrs = corrs.reset_index().rename({'index': 'Biomarker'}, axis=1)
    
    return corrs

def get_preselected_biomarkers(gene,
                               biomarkers,
                               sample_subset,
                               gene_effect_matrix, 
                               biomarker_matrix):
    corrs = quick_corr(
        gene_effect_matrix.loc[sample_subset, gene], biomarker_matrix.loc[:, biomarkers], 
        n=None
    )
    
    corrs = corrs.reset_index().rename({'index': 'Biomarker'}, axis=1)
    
    return corrs

def run_correlations_paired(pairs, mtx1, mtx2, method='Pearson'):
    results = []
    for p1, p2 in tqdm(pairs):
        aln_x, aln_y = mtx1[p1].dropna().align(mtx2[p2].dropna(), join='inner')
        n = len(aln_x)
        if aln_x.var() == 0 or aln_y.var() == 0 or n < 2:
            results.append(pd.Series({'Feature 1': p1, 'Feature 2': p2, f'{method} Correlation': np.nan, 'p-value': np.nan, 'n': n}))
        else:
            if method.lower() == 'pearson':
                corr, p = scipy.stats.pearsonr(aln_x, aln_y)
            elif method.lower() == 'spearman':
                corr, p = scipy.stats.spearmanr(aln_x, aln_y)
            results.append(pd.Series({'Feature 1': p1, 'Feature 2': p2, f'{method} Correlation': corr, 'p-value': p, 'n': n}))
    return pd.DataFrame(results)

# two-class comparison functions
def run_mann_whitney_u(full_matrix, group1, group2, group1_name=None, group2_name=None, **test_kws):
    if group1_name is None:
        group1_name = 'Group1'
    if group2_name is None:
        group2_name = 'Group2'
    
    mwu_results = scipy.stats.mannwhitneyu(full_matrix.reindex(index=group1), full_matrix.reindex(index=group2), **test_kws)
    
    output_df = pd.concat([
        (~full_matrix.reindex(index=group1).isna()).sum().reindex(index=full_matrix.columns).rename(f'{group1_name} n'),
        (~full_matrix.reindex(index=group2).isna()).sum().reindex(index=full_matrix.columns).rename(f'{group2_name} n'),
        pd.Series(mwu_results.statistic, index=full_matrix.columns).rename('U1'),
        pd.Series(mwu_results.pvalue, index=full_matrix.columns).rename('p')
    ], axis=1)
    
    fdr = pd.Series(scipy.stats.false_discovery_control(output_df.dropna()['p']), index=output_df.dropna().index).rename('FDR')
    output_df.loc[fdr.index, 'FDR'] = fdr
    
    return output_df

def run_ttest(full_matrix, group1, group2, group1_name=None, group2_name=None, **test_kws):
    if group1_name is None:
        group1_name = 'Group1'
    if group2_name is None:
        group2_name = 'Group2'
    
    t_results = scipy.stats.ttest_ind(full_matrix.reindex(index=group1), full_matrix.reindex(index=group2), **test_kws)
    
    output_df = pd.concat([
        full_matrix.reindex(index=group1).mean().rename(f'{group1_name} Mean'),
        full_matrix.reindex(index=group1).var().rename(f'{group1_name} Var'),
        (~full_matrix.reindex(index=group1).isna()).sum().rename(f'{group1_name} n'),
        full_matrix.reindex(index=group2).mean().rename(f'{group2_name} Mean'),
        full_matrix.reindex(index=group2).var().rename(f'{group2_name} Var'),
        (~full_matrix.reindex(index=group2).isna()).sum().rename(f'{group2_name} n'),
        (full_matrix.reindex(index=group1).mean() - full_matrix.reindex(index=group2).mean()).rename(f'Mean Difference'),
        pd.Series(t_results.statistic, index=full_matrix.columns).rename('t'),
        pd.Series(t_results.pvalue, index=full_matrix.columns).rename('p')
    ], axis=1)
    
    fdr = pd.Series(scipy.stats.false_discovery_control(output_df.dropna()['p']), index=output_df.dropna().index).rename('FDR')
    output_df.loc[fdr.index, 'FDR'] = fdr
    
    return output_df

def run_fishers_exact(full_matrix, group, outgroup):    
    test_inputs = pd.concat([
        full_matrix.reindex(index=group).sum().rename('k'), # observed successes
        (~full_matrix.reindex(index=group + outgroup).isna()).sum().rename('M'), # size of the population
        full_matrix.reindex(index=group + outgroup).sum().rename('n'), # possible successes in the population
        (~full_matrix.reindex(index=group).isna()).sum().rename('N'), # draws from the population
    ], axis=1)
    
    # survival function gives the probability of at least as many observed successes
    pvals = pd.Series(scipy.stats.hypergeom.sf(
        test_inputs["k"].astype(np.int64), 
        test_inputs["M"].astype(np.int64), 
        test_inputs["n"].astype(np.int64), 
        test_inputs["N"].astype(np.int64), 
        loc=1
    ), index=test_inputs.index).rename('p')
    
    fdrs = pd.Series(scipy.stats.false_discovery_control(pvals.dropna()), index=pvals.dropna().index).rename('FDR')
    
    test_inputs.loc[pvals.index, 'p'] = pvals
    test_inputs.loc[fdrs.index, 'FDR'] = fdrs
    
    return test_inputs

def run_mwu_and_ttest(full_matrix, group1, group2, group1_name=None, group2_name=None, equal_var=True, join_gene_metadata=None,
                      fishers_exact_matrix=None, **common_test_kws):
    mwu_res = run_mann_whitney_u(
        full_matrix, group1, group2, group1_name, group2_name, **common_test_kws
    )
    
    t_res = run_ttest(
        full_matrix, group1, group2, group1_name, group2_name, equal_var=equal_var, **common_test_kws
    )
    
    combined_res = t_res.rename({'p': 't p', 'FDR': 't FDR'}, axis=1).join(
        mwu_res.rename({'p': 'MWU p', 'FDR': 'MWU FDR'}, axis=1).loc[:, ['U1', 'MWU p', 'MWU FDR']]
    )
    
    if fishers_exact_matrix is not None:
        fe_res = run_fishers_exact(fishers_exact_matrix, group1, group2)
        fe_res[f'{group1_name} Dep'] = fe_res['k'] 
        fe_res[f'{group2_name} Dep'] = fe_res['n'] - fe_res['k']
        combined_res = combined_res.join(fe_res[[f'{group1_name} Dep', f'{group2_name} Dep', 'p', 'FDR']].rename({'p': 'Fisher p', 'FDR': 'Fisher FDR'}, axis=1))
    
    if join_gene_metadata is not None:
        combined_res = combined_res.join(join_gene_metadata)
    
    return combined_res

## gene set enrichment analyses
def merge_dicts(dict_list):
    final_dict = dict()
    for dct in dict_list:
        for k in dct:
            final_dict[k] = dct[k]
    return final_dict

def run_single_hypergeometric(genes_to_test, geneset_genes, all_genes, report_genes=False, 
                              preprocess_query=True, preprocess_geneset=True):
    original_geneset_size = len(geneset_genes)
    
    # preprocess sets to ensure validity in the universe set
    if preprocess_query:
        genes_to_test = [x for x in genes_to_test if x in all_genes]
    if preprocess_geneset:
        geneset_genes = [x for x in geneset_genes if x in all_genes]
    
    overlap = list(set(genes_to_test).intersection(set(geneset_genes)))
    
    success = len(overlap)
    total = len(all_genes)
    max_successes = len(geneset_genes)
    draws = len(genes_to_test)
    
    pval = scipy.stats.hypergeom.sf(success, total, max_successes, draws, loc=1)
    
    res = pd.Series({
        'n_overlap': success, 
        'effective_set_size': max_successes, 
        'original_set_size': original_geneset_size,
        'odds_ratio': ((success * (total - max_successes - draws + success)) / ((max_successes - success) * (draws - success))),
        'pval': pval})
    if report_genes:
        res.loc['overlap'] = ';'.join(overlap)
    return res

def run_multiple_hypergeometric(query_list, genesets, all_genes, report_genes=False):
    # preprocess query set
    query_clean = [x.split(" ")[0] for x in query_list]
    query_clean = [x for x in query_clean if x in all_genes]
    
    tests = []
    for gs in genesets:
        tests.append(
            run_single_hypergeometric(
                query_clean, genesets[gs], all_genes, report_genes=report_genes,
                preprocess_query=False, preprocess_geneset=True
            ).rename(gs)
        )
    tests = pd.concat(tests, axis=1).T
    tests['n_overlap'] = tests['n_overlap'].astype(int)
    tests['effective_set_size'] = tests['effective_set_size'].astype(int)
    tests['original_set_size'] = tests['original_set_size'].astype(int)
    tests['odds_ratio'] = tests['odds_ratio'].astype(float)
    tests['pval'] = tests['pval'].astype(float)
    tests['FDR'] = scipy.stats.false_discovery_control(tests['pval'])
    return tests.sort_values('FDR')

def run_rank_hypergeometric(ordered_series, genesets, 
                            n_top=1000, direction='left', 
                            all_genes=None, report_genes=False):
    if direction == 'left':
        query_genes = ordered_series.head(n_top).index.tolist()
    elif direction == 'right':
        query_genes = ordered_series.tail(n_top).index.tolist()
    elif direction == 'both':
        query_genes = ordered_series.head(n_top // 2).index.tolist() + ordered_series.tail(n_top // 2).index.tolist()
    else:
        raise ValueError(f'`direction` = "{direction}" not recognized. Use one of ["left", "right", "both"]')
    
    if all_genes is None:
        return run_multiple_hypergeometric(
            query_genes, genesets, 
            all_genes=ordered_series.index.str.split(' ').str[0].tolist(), 
            report_genes=report_genes
        )
    else:
        return run_multiple_hypergeometric(
            query_genes, genesets, 
            all_genes=all_genes, 
            report_genes=report_genes
        )

