#!/usr/bin/env python3

import os
import json
import subprocess
import pandas as pd
from Bio import SeqIO


def check_for_pseudogenes(df_blast, tID_col=None, fasta_file=None):
    """
    Identify Pseudo genes from Fasta Header
    param df_blast <pd.DataFrame>: df containing target IDs that should be checked
    param tID_col <str>: name of column in which tIDs are stored
    return <list>: is_pseudo_gene flags [True/False]
    """
    
    tIDs = df_blast[tID_col].unique()
    
    dict_pseudo = dict()
    for record in SeqIO.parse(fasta_file, 'fasta'):
        if record.id in tIDs:

            if record.id in dict_pseudo:
                print('tID {} already present'.format(record.id))

            if '[pseudo=true]' in record.description:
                dict_pseudo[record.id] = True
            else:
                dict_pseudo[record.id] = False
    
    col_pseudo = [dict_pseudo[tID] for tID in df_blast[tID_col]]

    return col_pseudo


def process_blast_result(blast_result_file, fasta_file, precomputed_table=None):
    """
    Import and process Blast result table
    :param df <pd.DataFrame>: Blast result table
    :param precomputed_table <str>: path to precomputed table; only implemented for demonstration purpose
    :return <pd.DataFrame>: processed Blast result table
    """

    if precomputed_table:
        if os.path.isfile(precomputed_table) == True:
            print('    Precomputed and processed BLAST result table found')
            print('        Importing table')
            return pd.read_table(precomputed_table, compression='gzip')
        else:
            print('    WARNING: precomputed blast table does not exist.')
            print('         Computing table...')

    not_precomputed = True
    df = pd.read_table(blast_result_file, header=None)
    columns = ['qseqid', 'sseqid', 'pident', 'nident', 'length', 'qstart', 'qend', 'evalue']
    df.columns = columns

    # calulate perc_cov for all queries
    col_coverage, col_length = calculate_perc_coverage(df, tID_col='qseqid', fasta_file=fasta_file)
    df['perc_coverage'] = col_coverage
    df['length'] = col_length

    # check for pseudo and add column
    df['pseudo'] = check_for_pseudogenes(df, tID_col='qseqid', fasta_file=fasta_file)
    
    if not_precomputed == True and precomputed_table:
        print('    Exporting precomputed and processed BLAST result table')
        df.to_csv(precomputed_table, sep='\t', index=False, compression='gzip')
        print('        Export done\n')

    return df


def calculate_perc_coverage(df_blast, tID_col=None, fasta_file=None):
    """
    Calculate Percent coverage
    :param df_blast <pd.DataFrame>:
    :param tID_col <str>: Name of target id columns
    :param fasta_file <str>: Path to fasta file
    :return <tuple>: coverages and sequence lengths
    """

    tIDs = df_blast[tID_col].unique()
    dict_length = _get_lengths(tIDs, fasta_file)

    col_length = list()
    col_coverage = list()

    for i, row in df_blast.iterrows():
        tID = row['qseqid']

        qstart = int(row['qstart'])-1
        qend = int(row['qend'])
        len_align = qend - qstart

        qseq_length = dict_length[tID]
        perc_coverage = len_align/ qend
        col_coverage.append(perc_coverage)
        col_length.append(qseq_length)
    
    return (col_coverage, col_length)


def _get_lengths(tIDs, fasta_file):
    """
    Get gene length from fasta file
    :param tIDs: list of target_ids found in fasta_file
    :param fasta_file: path to fasta file containing information
    :return: dict {tID: length}
    """

    dict_length = dict()
    for record in SeqIO.parse(fasta_file, 'fasta'):

        if record.id in tIDs:

            if record.id in dict_length:
                print('tID {} already present'.format(record.id))

            qseq_length = len(record.seq)
            dict_length[record.id] = qseq_length
    
    return dict_length


def gene_clusters_from_blast(df, sort_qseqs=False):
    """
    Create cluster genes based on Blast results table
    :param df <pd.DataFrame>: processed blast result table
    :return <dict>: gene_group to member assignment
    """

    df = df.copy()

    pseudo_genes = list(df[df['pseudo'] == True])

    # dict stores gene group identifiers and members
    dict_genegroups = dict()

    if sort_qseqs == True:
        qseqids = sort_qseqids_by_length(df)
    else:
        qseqids = list(df['qseqid'].unique())

    for qseqid in qseqids:
        # read blast filtered query
        list_sseqids = list()
        list_sseqids = df[df['qseqid'] == qseqid]['sseqid']

        # update list of seen genes
        clustered_genes = set([gene for group, genes in dict_genegroups.items() for gene in genes])
        seen_genes = [sseqid for sseqid in list_sseqids if sseqid in clustered_genes]

        # qseqs define a gene group {qseq: [associated sseqids]}
        # pseudo genes include truncated genes or gene with frame shifts and  should not define a gene cluster
        if qseqid not in pseudo_genes:

            # check if sseq are already part of an existing VF gene cluster
            if seen_genes:
                groupIDs = set([groupID for seen_gene in seen_genes \
                                for groupID, members in dict_genegroups.items() \
                                if seen_gene in members])

                # add sseqs to existing gene cluster
                groupID = sorted(list(groupIDs))[0]
                dict_genegroups[groupID].extend(list_sseqids)

            # all associated sseqs are no part of an existing cluster
            else:
                # create a new VF gene cluster
                if qseqid not in dict_genegroups.keys():
                    dict_genegroups[qseqid] = list()

                dict_genegroups[qseqid].extend(list_sseqids)

    return dict_genegroups


def generate_gene_cluster(dict_genegroups, faa_file, export_file=None):
    """
    Convert qseq_ids to cluster_ids and add gene description
    :param dict_genegroups <dict>: gene group id to member assignment
    :param export_file <str>: Path to json file
    :return <dict>: updated VF gene cluster db
    """

    tIDs = [tID for tID in dict_genegroups.keys()]
    qseqid2descr = _get_descriptions(tIDs, faa_file)

    counter = 1
    dict_clusters = dict()
    for qseqid, sseqids in dict_genegroups.items():
        cluster_id = 'strepVF_{0:0>4}'.format(counter)
        description = qseqid2descr[qseqid]
        dict_clusters[cluster_id] = {'members': sseqids, 'description':description}
        counter += 1

    if export_file:

        with open(export_file, 'w') as db:
            json.dump(dict_clusters, db)
            print('++ cluster db exported ++')

    return dict_clusters


def _get_descriptions(tIDs, fasta_file):
    """
    Extract description from fasta header
    :param tIDs <list>: list of target_ids found in fasta_file
    :param fasta_file <str>: Path to fasta file containing information
    :return <dict>: {tID: description}
    """

    dict_descr = dict()
    for record in SeqIO.parse(fasta_file, 'fasta'):

        if record.id in tIDs:

            if record.id in dict_descr:
                print('tID {} already present'.format(record.id))
            
            descr = record.description.split('[protein=')[1].split(']')[0]
            dict_descr[record.id] = descr
    
    return dict_descr


def sort_qseqids_by_length(df_blast):
    """
    Sort blast results by query seq gene length
    :param df_blast <list>: target IDs
    :return <list>: sorted query seq ids
    """

    # only unique qseqids
    mask = df_blast[['qseqid', 'length']].duplicated()
    # sort seqs by gene length
    sorted_qseqids = list(df_blast[~mask].sort_values('length', ascending=False)['qseqid'])

    return sorted_qseqids


def main():
    num_threads = 8

    base_dir = os.getcwd()
    faa_file = './data/fig4/streptococcus_VF.faa'
    blastdb_file = './data/BLAST/streptococcus_REF_PAN'
    blast_result_file = './data/fig4/BLASTP_VF_ident50.tsv.gz'
    blast_result_precomputed_file = './data/fig4/BLASTP_VF_ident50_precomputed.tsv.gz'

    # Virulence Factor gene cluster database
    db_file = './data/fig4/streptococcus_vir_db.json'

    blast_dir = '/usr/local/ncbi/blast/bin'
    blastp_bin = '{blast_dir}/blastp'.format(blast_dir=blast_dir)

    # if not already done run blastp virulence genes vs streptococcus database
    if os.path.isfile(blast_result_file) == False:
        listBlastParam = ['-evalue', '1E-6', '-num_threads', str(num_threads),
                      '-outfmt', '6 qseqid sseqid pident nident length qstart qend evalue']

        subprocess.check_call([blastp_bin, '-query', faa_file,
                               '-db', blastdb_file, '-out', blast_result_file] + listBlastParam)

    # NOTE: HERE function only imports a precomputed version of the processed BLAST result table (run time reduction)
    df_blast = process_blast_result(blast_result_file, faa_file, precomputed_table=blast_result_precomputed_file)

    # filter by percent coverage and percent identity
    print('ALL entries:', df_blast.shape[0])
    df_blast = df_blast[df_blast['perc_coverage'] >= 0.8]
    print('Filter coverage: ', df_blast.shape[0])
    df_blast = df_blast[df_blast['pident'] >= 50]
    print('Filter identity: ', df_blast.shape[0])
    
    #dict_genegroups = gene_clusters_from_blast(df_blast)
    dict_genegroups = gene_clusters_from_blast(df_blast)

    # remove duplicate seqs from gene clusters
    dict_genegroups = {group: list(set(members)) for group, members in dict_genegroups.items() }
    all_genes = len(df_blast['sseqid'].unique())
    groups = len(dict_genegroups)
    members = len(set([member for group, members in dict_genegroups.items() for member in members]))
    pseudos = len(df_blast[df_blast['pseudo'] == True])
    
    print('# of genes:\t', all_genes)
    print('# of pseudos:\t', pseudos)
    print('# of groups:\t', groups)
    print('# of members:\t', members)

    # export cluster to json file
    dict_cluster = generate_gene_cluster(dict_genegroups, faa_file, export_file=db_file)

if __name__ == '__main__':
    main()