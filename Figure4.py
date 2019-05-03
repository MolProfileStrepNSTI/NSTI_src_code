#!/usr/bin/env python3

import json
import os
import pandas as pd

from Bio import SeqIO
import VF_gene_clustering as IC
from DataOrganization import DataOrganization


def _get_descriptions(tIDs, fasta_file):
    """
    param tIDs: list of target_ids found in fasta_file
    param fasta_file: path to fasta file containing information
    returns: dict {tID: description}
    """

    dict_descr = dict()
    for record in SeqIO.parse(fasta_file, 'fasta'):
        if record.id in tIDs:
            if record.id in dict_descr:
                print('tID {} already present'.format(record.id))
            descr = record.description.split('[protein=')[1].split(']')[0]
            dict_descr[record.id] = descr
    return dict_descr


def assign_TPM2groups(df_VF, DO, genus=None):
    """
    table: geneGroup VFclass VFgene species_in_VFDB species_groupMembers [sample_TPMs]
    aggregated table
    :param df_VF <pd.DataFrame>:
    :param genus <str>: Genus Name
    """
    #dfVF [index(virCluster), 'VFclass', 'VF_gene_names', 'groupIDs']
    columns = ['VFclass', 'VFgene_names', 'groupIDs']
    df_VF = df_VF[columns]

    df_group2VFgene = dict()
    species = set([spec for k, v in DO.speciesVF_occurance.items() for spec in v])
    for i, row in df_VF.iterrows():
        VFclasses = row['VFclass']
        VFgene_names = row['VFgene_names']

        groupIDs = row['groupIDs'].split('|')
        if not groupIDs:
            print('no groupID found for', VFgene_names)
            continue

        #gID:  {VFclass: '', VFgene_names: '', species_in_VFDB: 0/1}
        for groupID in groupIDs:

            if groupID != '':

                if groupID not in df_group2VFgene:
                    df_group2VFgene[groupID] = {'VFclass': VFclasses, 'VFgene_names': VFgene_names}
                elif groupID not in df_group2VFgene and df_group2VFgene[groupID]['VFgene_names'] == VFgene_names:
                    continue
                else:
                    df_group2VFgene[groupID]['VFgene_names'] += '|{}'.format(VFgene_names)
                    df_group2VFgene[groupID]['VFclass'] += '|{}'.format(VFclasses)
                    print(groupID, VFclasses, VFgene_names, df_group2VFgene[groupID])
                
                # add binary if VFDB gene present in species
                for spec in species:
                    VF_in_spec = 0
                    for VFgene_name in VFgene_names.split('|'):
                        VF_in_spec += DO.speciesVF_occurance[VFgene_name][spec]
                        
                    df_group2VFgene[groupID]['{}_in_VFDB'.format(spec)] = VF_in_spec

    # add TPMs
    df_tpm = DO._import_tpm_db(to_df=True)
    info_cols = list(df_tpm.columns[:9])
    samples = list(df_tpm.columns[14:])
    
    df = pd.DataFrame.from_dict(df_group2VFgene, orient='index')
    vir_gIDs = list(df.index)
    species = list(species)

    # check if species genes in geneGroup
    for spec in species:
        col_name = '{}_groupMembers'.format(spec)
        df[col_name] = DO.identify_species_gene_clusters(vir_gIDs, species_name=spec, genus=genus)

    df = pd.concat([df, df_tpm.loc[vir_gIDs, info_cols+samples]], axis=1)
    df.sort_values('VFgene_names', inplace=True)

    return df


def create_VFfaa(tIDs, genus='streptococcus', faa_file=None):
    """
    Creates protein FASTA containing VirDB genes
    :param tIDs <list>: target Íds
    :param genus <str>: Genus Name
    :param faa_file <str>: Path to Virulence Factor fasta file
    """

    genus = genus.capitalize()
    ref_faa = '{data_dir}/{gen}_REF_PAN.faa'.format(data_dir=DO.data_dir, gen=genus.lower())

    if not faa_file:
        faa_file = '{data_dir}/{gen}_VF.faa'.format(data_dir=DO.data_dir, gen=genus.lower())

    faa_records = list()
    for record in SeqIO.parse(ref_faa, 'fasta'):

        if record.id in tIDs:
            faa_records.append(record)
    
    print('+++++ exporting protein fasta {faa}'.format(faa=faa_file))
    SeqIO.write(faa_records, faa_file, 'fasta')


def create_src_data_table(df_VFcounts, DO):
    """
    Create source data table
    :param df_VFcounts <pd.DataFrame>: 
    :param DO <DataOrganization object>: Initialized databases
    :return <pd.DataFrame>: Figure 4B source data table
    """

    # names in Fig4
    VF_targets = {'bca': ['bca'], 'cba': ['cba'], 'grab': ['grab'], 'sclA': ['sclA'],
                'sclB': ['sclB'], 'emm': ['emm', 'enn'], 'eno': ['eno'], 'lmb': ['lmb'], 
                'prtF1': ['prtF1', 'cpa'], 'prtF2': ['prtF2', 'fbaA'], 'sfbII': ['sfbII'], 'sfbX': ['sfbX'], 
                'fbp54': ['fbp54'], 'fbsB': ['fbsB'], 'sic': ['sic'], 'ska/skc/skg': ['ska'], 
                'endoS': ['endoS'], 'pilA': ['pilA'], 'scpA': ['scpA'], 
                'cepA (SPYCEP)': ['cepA'], 'IdeS': ['ideS'], 'slo': ['slo'],  
                'sagA': ['sagA'], 'spyA': ['spyA'], 'speA': ['speA'], 
                'speB': ['speB'], 'speC|speJ': ['speC', 'speJ'], 
                'speG': ['speG'], 'speH': ['speH'], 'speI': ['speI'], 
                'speK|speM': ['speK', 'speM'], 'smeZ': ['smeZ'], 'cfb': ['cfb']}
    
    df_tpm = DO._import_tpm_db(to_df=True)
    samples = ['2006', '2066', '3004', '3005', '3017', '3019', '3026', '3032', '3037', 
               '6007', '6010', '6017', '6023', '6025', '6026', '6033', '6057']

    temp = list()
    for name, VFs in VF_targets.items():
        x = df_VFcounts[df_VFcounts['VFgene_names'].str.contains('|'.join(VFs)) == True][samples].transpose().sum(axis=1)

            
        # manually added by locus Tag to genegroup assignment
        if name == 'cepA (SPYCEP)':
            group_id = 'STREP_1242'
            x = df_tpm.loc[group_id, samples]

        if name == 'prtF2':
            x = df_VFcounts[df_VFcounts['VFgene_names'].str.contains('|'.join(VFs)) == True][samples]
            # manually curated, belongs to cba
            x = x.drop('STREP_5086')
            x = x.transpose().sum(axis=1)
        
        x.name = name
        temp.append(x)

    df = pd.DataFrame(temp).transpose()
    df.index = df.index.astype(int)
    return df


if __name__ == '__main__':
    data_dir = 'data'
    result_dir = 'results/figure4'
    os.makedirs(result_dir,exist_ok=True)

    genus = 'Streptococcus'
    reanalyze = True
    create_new_fasta = False

    print('***** {}'.format(genus))
    DO = DataOrganization(data_dir)
    
    # get all pyogenes and agalactiae genes from VFDB and create list of associated tIDs
    # create FASTA file as input for BLAST-based geneclustering
    faa_file = '{data_dir}/{gen}_VF.faa'.format(data_dir=DO.data_dir, gen=genus.lower())

    if os.path.isfile(faa_file) == False or reanalyze == True:
        print('\n+++++ STARTING BLAST Virulence Clustering')

        if create_new_fasta == True:
            print('\n+++++ Creating FASTA')
            DO._import_VFDB(genus, species_filter=['pyogenes', 'agalactiae'])
            geneinfo_db = DO._import_geneinfo_db(genus)
            tIDs = list()

            for VFgene_name, locTs in DO.dict_vir.items():
                #convert locTs to tID
                tIDs.extend([tID for locT in locTs for tID in geneinfo_db if locT == geneinfo_db[tID]['locT']])

            tIDs = list(set(tIDs))
            create_VFfaa(tIDs, genus='streptococcus', faa_file=faa_file)

        # run blastp clustering
        print('+++++ RUN BLAST')
        IC.main()
        print('+++++ CLUSTERING done')

    DO.import_databases(genus, species_filter=['pyogenes', 'agalactiae'])
    print('+++++ fetching VF cluster associations +++++')
    # add VF cluster description and VF statistics
    DO.get_VFDB_names(genus)

    # create and export VF cluster data table
    dfVF_cluster = pd.DataFrame.from_dict(DO.VFcluster_db, orient='index')
    dfVF_cluster.sort_values('VFgene_names', inplace=True)
    #vir_expr_file = 'data/fig4/{gen}_VFcluster_data_+manual_dysgenes.xlsx'.format(gen=genus)
    vir_expr_file = f'{result_dir}/Fig4B_src.xlsx'
    writer = pd.ExcelWriter(vir_expr_file, engine='xlsxwriter')
    #dfVF_cluster.to_excel(writer, sheet_name='VF_cluster')

    print('+++++ assigning TPMs +++++')
    df_VFcounts = assign_TPM2groups(dfVF_cluster, DO, genus=genus)
    df_src = create_src_data_table(df_VFcounts, DO)
    df_src.to_excel(writer, sheet_name='Fig4B_src')

    writer.save()