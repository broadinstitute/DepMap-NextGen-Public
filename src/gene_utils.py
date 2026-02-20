import pandas as pd
import numpy as np
import scipy
from tqdm import tqdm
from data_utils import *

hgnc_annotations = load_file('hgnc_table.csv')
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

# helper function for resolving multimapping annotations
def combine_annotations(binary_matrix, annotation_columns=[], use_column_name=True):
    annotations = pd.Series(np.nan, index=binary_matrix.index)
    for c in annotation_columns:
        to_replace = list(set(binary_matrix[c].loc[lambda x: x].index.tolist()) & # is annotated in the next column
                          set(annotations.isna().loc[lambda x: x].index.tolist())) # has not already been annotated
        if use_column_name:
            annotations.loc[to_replace] = c
        else:
            annotations.loc[to_replace] = binary_matrix.loc[to_replace, c]
    return annotations

def expand_model_matrix_to_screens(matrix, model_mapping):
    return matrix.merge(
        model_mapping.to_frame(), how='inner', left_index=True, right_index=True
    ).set_index(model_mapping.name, drop=True)

def merge_dicts(dict_list):
    final_dict = dict()
    for dct in dict_list:
        for k in dct:
            final_dict[k] = dct[k]
    return final_dict

def geneset_min_coverage_similarity(genesets_to_genes_dict, geneset_subset=None):
    if geneset_subset is None:
        geneset_subset = genesets_to_genes_dict.keys()
    min_coverage_similarity = np.zeros((len(geneset_subset), len(geneset_subset)))
    
    for i, gsi in enumerate(geneset_subset):
        for j, gsj in enumerate(geneset_subset):
            intersection = set(genesets_to_genes_dict[gsi]) & set(genesets_to_genes_dict[gsj])
            min_size = min(len(genesets_to_genes_dict[gsi]), len(genesets_to_genes_dict[gsj]))
            min_coverage_similarity[i, j] = len(intersection) / min_size
    
    return pd.DataFrame(min_coverage_similarity, index=geneset_subset, columns=geneset_subset)

def connected_component(adjacency_matrix, source):
    if type(adjacency_matrix) == pd.DataFrame:
        assert (source in adjacency_matrix.columns.tolist())
        adj_mtx = adjacency_matrix.values
        src = adjacency_matrix.columns.tolist().index(source)
    elif type(adjacency_matrix) == np.ndarray:
        assert (source in range(adjacency_matrix.shape[0]))
        adj_mtx = adjacency_matrix
        src = source
    
    seen = set()
    frontier = set([src])
    while len(frontier) > 0:
        visited = frontier.pop()
        seen.add(visited)
        neighbors = np.where(adj_mtx[:, visited] == 1)[0].tolist()
        frontier = frontier.union(set(neighbors) - seen)
    if type(adjacency_matrix) == pd.DataFrame:
        return sorted([adjacency_matrix.columns.tolist()[i] for i in seen])
    else:
        return sorted(list(seen))
    
def partition_by_connected_component(adjacency_matrix):
    if type(adjacency_matrix) == pd.DataFrame:
        all_columns = set(adjacency_matrix.columns.tolist())
    elif type(adjacency_matrix) == np.ndarray:
        all_columns = set(range(adjacency_matrix.shape[0]))
    partitions = []
    while len(all_columns) > 0:
        visited = all_columns.pop()
        cc = connected_component(adjacency_matrix, visited)
        partitions.append(cc)
        all_columns = all_columns - set(cc)
    clusters = dict()
    for i, grp in enumerate(partitions):
        for el in grp:
            clusters[el] = i + 1
    return pd.Series(clusters)

def clean_geneset_name(string, nth=2, remove_prefix=True):
    string = string.replace('_', ' ')
    split_str = string.split(' ')[remove_prefix:]
    new_str_list = []
    for j, x in enumerate(split_str):
        if j % nth == 0:
            new_str_list.append('\n')
        else:
            new_str_list.append(' ')
        new_str_list.append(x)
    return ''.join(new_str_list).strip()
