#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os

import networkx as nx
import numpy as np
import pandas as pd
from scipy.stats import mstats
from statsmodels.stats.multitest import fdrcorrection

import obo
from DataOrganization import DataOrganization


def create_gml(O, G, GML_file='results/figure3/induced_network.gml', lst_strep=None,  
               lst_poly=None, lst_shared=None):
    """
    Create and export induced GO network in gml file format
    :param O <object>: GO ontology NetworkX graph
    :param G <object>: GO ontology NetworkX subgraph
    :param GML_file <str>: Path to output gml file
    :param lst_strep <list>: induced GO terms in Streptococcus group
    :param lst_poly <list>: induced GO terms in Polymicrobial group
    :param lst_shared <list>: GO terms shared between both groups
    """
    
    results_dir, fname = GML_file.rsplit('/', 1)
    os.makedirs(results_dir, exist_ok=True)

    color = {'strep': '#FF0000',
              'poly': '#0000FF',
              'shared': '#E0E0E0',
              'edge': '#000000',
              'other': '#E0E0E0'}
    
    shape = {'strep': '"ellipse"',
              'poly': '"ellipse"',
              'shared': '"ellipse"',
              'other': '"circle"'}
    
    width = {'strep': '4',
          'poly': '4',
          'shared': '2',
          'other': '1'}
    
    with open(GML_file, 'w') as f:
        ### create graph
        f.write('graph\n[\n')
        f.write('\tdirected 1\n')
        f.write('\tIsPlanar 1\n')
        
        ### write nodes
        for n in G.nodes():
             
            if n in lst_strep:
                classification = 'strep'
            elif n in lst_poly:
                classification = 'poly'
            elif n in lst_shared:
                classification = 'shared'
            else:
                classification = 'other'
            
            ### open node
            f.write('\tnode\n')
            f.write('\t[\n')
            
            ### node info
            f.write('\t\tid {}\n'.format(int(n.split(':')[1])))
            f.write('\t\tname "{}"\n'.format(n))
            f.write('\t\tlabel "{}"\n'.format(O.node[n]['name']))
            f.write('\t\tclassification "{}"\n'.format(classification))
            
            ### graphics
            f.write('\t\tgraphics\n')
            f.write('\t\t[\n')            
            
            f.write('\t\t\twidth {}\n'.format(width[classification]))
            f.write('\t\t\ttype {}\n'.format(shape[classification]))
            f.write('\t\t\tfill "{}"\n'.format(color[classification]))
            f.write('\t\t\toutline "{}"\n'.format(color['edge']))
            f.write('\t\t]\n')
            
            ### close node
            f.write('\t]\n')
        
        ### write edges
        for e in G.edges():
            ###open edge
            f.write('\tedge\n')
            f.write('\t[\n')
            
            ### edge info
            f.write('\t\tsource {}\n'.format(int(e[0].split(':')[1])))
            f.write('\t\ttarget {}\n'.format(int(e[1].split(':')[1])))
            
            ### graphics
            f.write('\t\tgraphics\n')
            f.write('\t\t[\n')
            f.write('\t\t\twidth 1\n')
            f.write('\t\t\ttype "line"\n')
            f.write('\t\t\tfill "{}"\n'.format(color['edge']))
            f.write('\t\t]\n')
            
            ### close edge
            f.write('\t]\n')
        ### close graph
        f.write(']')


def distances_to_node(O, GO, target_node='GO:0008150', get_path=False):
    """
    Find longest and shortest path from GO id to target node 
    :param O <object>: GO ontology NetworkX graph
    :param GO <str>: start GO id
    :param tarhet_node <str>: target node
    :param get_path <bool>: Specify to return path (True) or distance (False)
    :return <tuple>: (shortest dist/path, longest dist/path) 
    """
    
    len_min = 255
    len_max = 0
    path_minlength = list()
    path_maxlength = list()
    
    for path in nx.all_simple_paths(O, GO, target_node):
        new_path_length = len(path)
        if new_path_length < len_min:
            len_min = new_path_length
            path_minlength = path
        if new_path_length > len_max:
            len_max = new_path_length
            path_maxlength = path

    if get_path == True:
        return path_minlength, path_maxlength
    else:
        return len_min, len_max


def create_subgraph(O, terms, length='min', target_node='GO:0008150'):
    '''
    :param O <object>: GO ontology NetworkX graph
    :param terms <list>: GO terms to extract
    :param length <str>: min/max build graph only with shortest/ longest paths to target_node
    :param target_node <str>: root node id of the subgraph
    :return: subgraph build with supplied GO terms
    '''

    nodes = set()
    edges = set()
    for term in terms:
        if term in nodes:
            continue
        short_path, long_path = distances_to_node(O, term, target_node=target_node, get_path=True)
        
        if length == 'min':
            nodes.update(short_path)
            short_edges = [(short_path[i], short_path[i+1]) for i in range(len(short_path)) if i < len(short_path)-1]
            edges.update(short_edges)
        if length == 'max':
            nodes.update(long_path)
            long_edges = [(long_path[i], long_path[i+1]) for i in range(len(long_path)) if i < len(long_path)-1]
            edges.update(long_edges)
    
    G = nx.DiGraph()
    G.add_nodes_from(list(nodes))
    G.add_edges_from(list(edges))
    return G

def compare_groups(df, group0, group1, group0_name='group0', group1_name='group1'):
    """
    Calculate log2 fold change and tests for statistical difference 
    (Kruskal-Wallis + FDR correction)
    :param <pd.DataFrame>: Table with normalized read counts for each GO_ID and sample
    :param group0 <list>: Samples of group 0
    :param group1 <list>: Samples of group 1
    :return <pd.DataFrame>: Table with added statistics
    """

    dict_results = {GOID: {'pval': 0.0, 'log2fc': 0.0, 'mean_{}_tpm'.format(group0_name): 0.0, 
                           'mean_{}_tpm'.format(group1_name): 0.0, 'padj': 0.0} for GOID in df.index}

    for GOID in df.index:
        GO_group0 = np.array(df.loc[GOID, group0])
        GO_group1 = np.array(df.loc[GOID, group1])
        
        try:
            H, pval = mstats.kruskalwallis(GO_group0, GO_group1)
        except:
            pval = np.nan
        dict_results[GOID]['pval'] = pval
        
        mean_group0 = np.mean(GO_group0)
        mean_group1 = np.mean(GO_group1)
        dict_results[GOID]['mean_{}_tpm'.format(group0_name)] = mean_group0
        dict_results[GOID]['mean_{}_tpm'.format(group1_name)] = mean_group1
        
        try:
            log2fc = np.log2(mean_group0/ mean_group1)
        except:
            log2fc = np.nan
        dict_results[GOID]['log2fc'] = log2fc
        
    df_results = pd.DataFrame.from_dict(dict_results, orient='index')
    df_results = df_results.replace([np.inf, -np.inf], np.nan)
    df_results = df_results.fillna(0.0)

    df_results = do_FDR_correction(df_results)

    return df_results


def do_FDR_correction(df):
    """
    Do FDR correction and add results to dataframe
    # code from gonenrich module of goenrich package by Jan Rudolph (jdrudolph)
    # https://github.com/jdrudolph/goenrich/blob/master/goenrich/enrich.py
    :param df <pd.DataFrame>: GO expression data
    :return <pd.DataFrame>: GO expression data w/ FDR results
    """
    
    _p = np.array(df['pval'])
    # create array of len corresponding to p
    padj = _p.copy()
    rej = _p.copy()
    # list of bools not nan
    mask = ~np.isnan(_p)
    # remove false entries
    p = _p[mask]
    _rej, _padj = fdrcorrection(p)
    # only change values not nan values
    rej[mask] = _rej
    padj[mask] = _padj

    df['padj'] = padj
    df['rejected'] = rej

    return df


def get_GO_descriptions(GO_ids, GO_db):
    """
    Retrieve GO names and Descriptions for specified GO ids
    :param GO_ids <list>: GO ids
    :param GO_db <dict>: DB w/ GO names and descriptions
    :return <tuple>: GO aspects, GO names
    """

    col_aspect = list()
    col_name = list()

    for GO in GO_ids:
        if GO in GO_db:
            aspect = GO_db[GO]['aspect']
            name = GO_db[GO]['name']
        else:
            aspect = ''
            name = ''
        
        col_aspect.append(aspect)
        col_name.append(name)
    
    return col_aspect, col_name


if __name__ == '__main__':
    results_dir = 'results/figure3'
    os.makedirs(results_dir, exist_ok=True)

    DO = DataOrganization('data')
    DO._import_sample_classification()
    dict_classification = DO.class_db

    dict_GOterm = DO.get_GO_terms()

    # update obsolete GO ids to recent versions
    dict_GOobsolete = {'GO:0004091': 'GO:0052689',
                       'GO:0006184': None,
                       'GO:0006467': 'GO:0003756',
                       'GO:0006200': None,
                       'GO:0007047': 'GO:0071555',
                       'GO:0008969': 'GO:0101006',
                       'GO:0009021': 'GO:0030697',
                       'GO:0009296': None,
                       'GO:0030092': None,
                       'GO:0043064': None,
                       'GO:0015300': 'GO:0015297',
                       'GO:0015563': 'GO:0022857',
                       'GO:0016023': 'GO:0031410',
                       'GO:0016044': 'GO:0061024',
                       'GO:0016820': 'GO:0042626',
                       'GO:0016876': 'GO:0004812',
                       'GO:0022891': 'GO:0022857',
                       'GO:0042963': 'GO:0019068',
                       'GO:0050827': 'GO:0090729',
                       'GO:0003840': 'GO:0036374'}

    # initialize GO ontology
    obo_file = DO.get_ontology()
    O = obo.ontology(obo_file)

    # initialize sample tpm table
    print('+++++ IMPORTING RAW DATA +++++')
    print()

    df_tpm = pd.read_excel('data/GO_normalized_reads.xlsx', index_col=0)
    # aggregate sample reads by GO id
    df_tpm = df_tpm.groupby('GO_Id').sum()
    samples = df_tpm.columns

    # only classified samples
    samples = [sample for sample in samples for cl_sample in dict_classification if cl_sample in sample]
    strep_samples = [sample for sample in samples if dict_classification[sample] == 'Streptococcus']
    poly_samples = [sample for sample in samples if dict_classification[sample] == 'Multiinfection']

    ### analysis
    threshold = 2.0
    min_tpm = 100
    filter_aspect = 'biological_process'

    print('+++++ Comparing groups +++++')
    df_compare = compare_groups(df_tpm, strep_samples, poly_samples, 
                                group0_name='strep', group1_name='poly')

    ### add GO aspect and term
    col_aspect, col_name = get_GO_descriptions(df_compare.index, dict_GOterm)
    df_compare['aspect'] = col_aspect
    df_compare['name'] = col_name

    columns = ['aspect', 'name', 'log2fc', 'mean_strep_tpm', 'mean_poly_tpm', 'pval', 
            'padj', 'rejected']
    df_compare = df_compare[columns]

    df_final = pd.concat([df_compare, df_tpm], join='inner', axis=1)
    columns = ['name', 'log2fc', 'mean_strep_tpm', 'mean_poly_tpm', 'padj']
    columns.extend(strep_samples)
    columns.extend(poly_samples)

    ## only BP
    df_compare = df_compare[(df_compare['aspect'] == filter_aspect)]
    df_final = df_final[(df_final['aspect'] == filter_aspect)]

    ## remove obsolete GOIDs
    lst_annotated = [GO for GO in df_compare.index if GO not in dict_GOobsolete]
    df_compare = df_compare.loc[lst_annotated, ]

    ## only consider GOs with a minimum TPM in any group
    df_compare = df_compare[(df_compare['mean_strep_tpm'] >= min_tpm) | (df_compare['mean_poly_tpm'] >= min_tpm)]
    lst_expressed = list(df_compare.index)

    ## only significantly enriched GOs
    lst_strep = list(df_compare[(df_compare['rejected'] == 1.0) & (df_compare['log2fc'] >= threshold)].index)
    lst_poly = list(df_compare[(df_compare['rejected'] == 1.0) & (df_compare['log2fc'] <= -threshold)].index)

    lst_remaining = [GO for GO in df_compare.index if GO not in lst_strep]
    lst_remaining = [GO for GO in lst_remaining if GO not in lst_poly]

    print()
    print('    Filtered aspect:', filter_aspect)
    print('    Filtered min. TPM:', min_tpm)
    print('    BP:', len(lst_annotated))
    print('    expressed GO terms:', len(lst_expressed))
    print('    Upredulated in Strep', len(lst_strep))
    print('    Upregulated in Poly', len(lst_poly))
    print('    unchanged GO terms', len(lst_remaining))
    print()

    """
    aspect: biological_process
    min_TPM: 100
    BP: 683
    expressed: 374
    strep 50
    poly 32
    remaining 292
    """

    induced_combined = lst_strep + lst_poly

    ### export src data for figure 3 (GO ID order corresponds to Figure 3)
    print('+++++ Exporting Source Data +++++')
    print()
    order_file = 'data/fig3/GO_id_order.txt'
    GO_order = list(pd.read_table(order_file, index_col=0, header=None).index)
    
    writer = pd.ExcelWriter('{rdir}/Fig3_src.xlsx'.format(rdir=results_dir), engine='xlsxwriter')
    # Table corresponds to Supplementary Table 7
    df_final.loc[lst_expressed, columns].to_excel(writer, sheet_name='Figure 3 src data')
    # Table corresponds to Figure 3 data order
    df_final.loc[GO_order, columns].to_excel(writer, sheet_name='Figure 3 data sorted')
    writer.save()

    ### export GO subnetwork
    print('+++++ Exporting induced GO subgraph +++++')
    G = create_subgraph(O, induced_combined, length='max', target_node='GO:0008150')

    create_gml(O, G, lst_strep=lst_strep, lst_poly=lst_poly, 
               lst_shared=lst_remaining, 
               GML_file='results/figure3/induced_network.gml')
    print('+++++ DONE +++++')
    print()