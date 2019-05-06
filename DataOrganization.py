#!/usr/bin/env python3

import os
import gzip
import json
import pandas as pd
import urllib.request

class DataOrganization:
    """
    Imports databases, sustaining INFECT compatibility
    """

    def __init__(self, data_dir):
        
        self.data_dir = data_dir

        self.assembly_db_file = '{data_dir}/Assembly_DB.json.gz'.format(data_dir=self.data_dir)
        self.groupmember_db_file = '{data_dir}/groupmember.json.gz'.format(data_dir=self.data_dir)
        
        self.tpm_db_file = '{data_dir}/Streptococcus_tpm_db.json.gz'.format(data_dir=self.data_dir)
        self.VFcluster_db_file = '{data_dir}/fig4/streptococcus_vir_db.json'.format(data_dir=self.data_dir)

        self.obo_file = '{data_dir}/go-basic.obo.gz'.format(data_dir=self.data_dir)
        self.class_db_file = '{data_dir}/sample_class_db.json.gz'.format(data_dir=self.data_dir)
        self.goterm_db_file = '{data_dir}/goterm_db.json.gz'.format(data_dir=self.data_dir)

        self.assembly_db = None
        self.geneinfo_db = None
        self.tpm_db = None
        self.class_db = None
        self.goterm_db = None


    def import_databases(self, genus, species_filter=None):
        """
        Imports VFDB for genus with compatible locT
        Imports VFcluster_db {VF_clusterid: {'members': [members], 'description': 'description'}}
        Imports groupmember_db {gene_group: [group_members]}
        :param genus <str>: Genus name
        :param species_filter <str>: Species name for data cleaning
        """

        self._import_VFDB(genus, species_filter=species_filter)
        self._import_VFcluster_db()
        self._import_groupmember_assignment(reverse=True)
        self._import_tpm_db()


    def _import_groupmember_assignment(self, reverse=False, return_db=False):
        """
        Import gene group to member assignment
        :param reverse <bool>: group to member or member to group assignment
        :param return_db <bool>: return database as dict
        returns: dict with gene to gene group assignment
        """

        with gzip.GzipFile(self.groupmember_db_file, 'r') as db:
            groupmember_db = json.loads(db.read())
        
        if reverse == True:
            groupmember_db = {member: group for group, members in groupmember_db.items() \
                                            for member in members}
        if return_db == True:
            return groupmember_db
        else:
            self.groupmember_db = groupmember_db


    def _import_locT_converter(self, genus):
        """
        :param genus <str>: genus name (not case sensitive) eg: Streptococcus
        return <dict>: locus Tag converter {RefSeq_locT: old_locT}
        """

        locT_converter_file = '{locT_dir}/{gen}_locTDatabase.json.gz'.format(locT_dir=self.data_dir, gen=genus.lower())
        
        with gzip.GzipFile(locT_converter_file, 'r') as db:
            locTconverter = json.load(db)

        return locTconverter


    def _import_tpm_db(self, to_df=False):
        """
        Import Expression data
        """

        if not self.tpm_db:
            with gzip.GzipFile(self.tpm_db_file, 'r') as db:
                self.tpm_db = json.loads(db.read())
        
        if to_df == True:
            df = pd.DataFrame.from_dict(self.tpm_db, orient='index')
            df.columns = [str(col) for col in df.columns]
            return df


    def _import_geneinfo_db(self, genus):
        """
        Import geneinfo database
        param genus: genus name (not case sensitive) eg: Streptococcus
        return <dict>: gene info database {tID: {'ascN': '', 'genus': '', 'locT': '', 'protein_id': '', 'description': ''}}
        """
        
        if not self.geneinfo_db:
            geneinfo_file = '{data_dir}/{gen}_GeneInfoDatabase.json.gz'.format(data_dir= self.data_dir, gen=genus.lower())
            
            with gzip.GzipFile(geneinfo_file, 'r') as db:
                self.geneinfo_db = json.load(db)
        
        return self.geneinfo_db


    def _import_refseqAssemblyDb(self):
        """
        Import RefSeq assembly db into class dict
        :return <dict>: RefSeq assembly db {AccN: 'Assembly_name': '', 'Genus': ...}
        """

        if not self.assembly_db:            
            with gzip.GzipFile(self.assembly_db_file, 'r') as db:
                self.assembly_db = json.loads(db.read().decode('UTF-8'))
        
        return self.assembly_db


    def _import_VFcluster_db(self):
        """
        Import database containing gene cluster to Virulence factor assignment
        """

        with open(self.VFcluster_db_file, 'r') as db:
            self.VFcluster_db = json.load(db)


    def _import_VFDB(self, genus, species_filter=None):
        """
        Import VFDB from xlsx file, converts locTs (compatible to RefSeq)
        
        :param genus <str>: Genus name for data cleaning
        :param species_filter <str>: Species name for data cleaning
        creates: dict_vir = {gene_name: [compatible_locTs]}
        creates: dict_vir2class = {gene_name: vir_class}
        """

        print('+++++ importing VFDB +++++')
        genus = genus.capitalize()
        VF_file = '{data_dir}/fig4/{gen}_VFs_comparsion.xls'.format(data_dir=self.data_dir, gen=genus.capitalize())
        #reverse locTconverter
        self.locTconverter = {j: i for i, j in self._import_locT_converter(genus).items()}
        
        ### exclude first line: contains Table name
        ### exclude last line: contains only version info
        df_vir = pd.read_excel(VF_file, header=1, skipfooter=1)
        df_vir.fillna('', inplace=True)
        
        if genus == 'Streptococcus':
            data_cleaning = {'R6 suface protein': 'R6_sfp','Capsule': 'has/neu/cps'}

            if species_filter == 'S. agalactiae':
                data_cleaning['Capsule'] = 'neu/cps'
                data_cleaning['Pilus island 1'] = 'PI-1'
                data_cleaning['Agglutinin receptor'] = 'AgI/II,ssp-5'
            else:
                data_cleaning['Agglutinin receptor'] = 'AgI/II,ssp-5'
                data_cleaning['Capsule'] = 'neu/cps/has'
                data_cleaning['Pilus island 1'] = 'PI-1'


        elif genus == 'Staphylococcus':
            data_cleaning = {'Capsule': 'cap'}
            #code = 'STAPH'
            if not species_filter:
                species_filter = 'S. aureus'

        elif genus == 'Escherichia':
            data_cleaning = {'LEE locus encoded TTSS': 'LEE',
                             'ACE T6SS': 'ACE',
                             'SCI-I T6SS': 'SCI-I'}
            #code = 'ESCHE'
        elif genus == 'Vibrio':
            data_cleaning = {'TTSS-1 secreted effectors': 'VOP',
                             'TTSS-2': 'TTSS-2'}
            #code = 'VIBRI'
        
        ### select included species
        species = df_vir.columns[2:]
        genomes = [info.split()[1] for info in df_vir.loc[0, species]]
        
        if species_filter:
            if isinstance(species_filter, list):
                selected_genomes = [gen for gen, spec in zip(genomes, species) for spec_filter in species_filter if spec_filter in spec]
                species_selected = [spec for spec in species for spec_filter in species_filter if spec_filter in spec]
            else:
                selected_genomes = [gen for gen, spec in zip(genomes, species) if species_filter in spec]
                species_selected = [spec for spec in species if species_filter in spec]
            print(species_selected)
        else:
            selected_genomes = [gen for gen, spec in zip(genomes, species)]
            species_selected = species
        
        if genus == 'Vibrio':
            species = df_vir.columns[2:]
            genomes = [info.split()[2] for info in df_vir.loc[0, species]]
            selected_genomes = [gen for gen, spec in zip(genomes, species)]
        
        df_selected = df_vir[(df_vir.loc[:, species_selected] != '').any(axis=1)]
        self.dict_vir = dict()
        
        virulence_factor = ''
        self.vir2class = dict()

        col_related_genes = list()
        for virulence_factor, gene_name in zip(df_selected['Virulence factors'], df_selected['Related genes']):
            if virulence_factor != '':
                VF = virulence_factor
            if gene_name != '':
                if gene_name == '-':
                    gene_name = data_cleaning[VF]
            col_related_genes.append(gene_name)
        df_selected.loc[:, 'Related genes'] = col_related_genes

        # create dict for species-specific virulence gene absence/ presence
        self.speciesVF_occurance = {vir_gene: {} for vir_gene in df_selected['Related genes']}
        for spec in species_filter:
            for vir_gene, present in zip(df_selected['Related genes'], df_selected[[spec_col for spec_col in species_selected if spec in spec_col]].any(axis=1)):
                self.speciesVF_occurance[vir_gene][spec] = int(present)

        for i, row in df_selected.iterrows():
            if row['Virulence factors'] != '':
                virulence_factor = row['Virulence factors']
            gene_name = row['Related genes']
            #if gene_name == '-':     #NOTE: related genes are already cleaned
            #    gene_name = data_cleaning[virulence_factor]
            if gene_name != '':
                self.vir2class[gene_name] = virulence_factor
                ascN2locT = {ascN: df_selected.loc[i, spec].replace(',', '').replace('*', '').split()\
                                            for spec, ascN in zip(species, genomes)}
                # get list of old locTs
                old_locTs = [locT for ascN, locTs in ascN2locT.items() for locT in locTs]
                new_locTs = [self.locTconverter[locT] for locT in old_locTs if locT in self.locTconverter]
                #creates dict(gene_name: [compatible_locTs])
                self.dict_vir[gene_name] = new_locTs
    

    def _import_sample_classification(self):
        """
        Import Sample classification
        """
        
        if not self.class_db:            
            with gzip.GzipFile(self.class_db_file, 'r') as db:
                self.class_db = json.loads(db.read().decode('UTF-8'))


    def get_VFDB_names(self, genus):
        """
        Assign VFDB names to VF cluster
        iterate VF cluster and assign corresponding checked VF genes to gene cluster associations
        :param genus <str>: Genus Name
        """
        
        for cluster_id in self.VFcluster_db:
            members = self.VFcluster_db[cluster_id]['members']
            VFclasses, VFgene_names  = self._get_virclasses(members, genus)
            groupIDs = set([self.groupmember_db[tID] for tID in members if tID in self.groupmember_db])
            self.VFcluster_db[cluster_id]['#members'] = len(members)
            self.VFcluster_db[cluster_id]['VFclass'] = '|'.join(sorted(list(VFclasses)))
            self.VFcluster_db[cluster_id]['VFgene_names'] = '|'.join(sorted(list(VFgene_names)))
            self.VFcluster_db[cluster_id]['groupIDs'] = '|'.join(sorted(list(groupIDs)))

            #print(cluster_id, list(groupIDs), len(members), VFclasses, VFgene_names, self.VFcluster_db[cluster_id]['description'])


    def _get_virclasses(self, members, genus):
        """
        Identify associated virulence classes/genes that include INFECT genes
        param members <list>: target IDs
        returns: list of VFclasses
        """
        
        geneinfo_db = self._import_geneinfo_db(genus)

        member_locTs = [geneinfo_db[tID]['locT'] for tID in members]
        #check if locT part of a VFclass
        VFgene_names = set([VFgene_name for VFgene_name, VF_locTs in self.dict_vir.items() 
                        for locT in member_locTs if locT in VF_locTs])
        VFclasses = set([self.vir2class[VFgene_name] for VFgene_name in VFgene_names])

        return (list(VFclasses), list(VFgene_names))


    def identify_species_gene_clusters(self, gIDs, genus=None, species_name='pyogenes'):
        """
        Identify Strep gene clusters containing target_species genes
        :param gIDs <list>: group IDs
        :param genus <str>: Genus name
        :param species_name <str>: specify to find species specific gene clusters
        """

        ## initialize variables
        if not genus:
            genus = species_name.split(' ')[0].capitalize()
            if genus == 'S.':
                genus = 'Streptococcus'
        gIDs = list(set(gIDs))
        #spec_cluster = {gID: False for gID in gIDs}
        spec_cluster = list()
        
        ## import databases
        Group2MemberConverter = self._import_groupmember_assignment(reverse=False, return_db=True)
        geneinfo_db = self._import_geneinfo_db(genus)
        assembly_db = self._import_refseqAssemblyDb()
        
        ## get all tIDs + assemblyIDs of gID members
        for gID in gIDs:
            gID_in_spec = 0

            if gID != '':
                group_members = Group2MemberConverter[gID]

                for member in group_members:
                    ## get GCF
                    ascN = geneinfo_db[member]['ascN']

                    for GCF in assembly_db:
                        ascNs = assembly_db[GCF]['accN']

                        if ascN in ascNs:
                            ## get organism name
                            organism = assembly_db[GCF]['organism']

                            if species_name in organism.lower():
                                #spec_cluster[gID] = True
                                gID_in_spec = 1
            
            spec_cluster.append(gID_in_spec)
    
        return spec_cluster


    def _VirlocTs2gIDs(self):
        """
        Virulence gene Locus Tags to GroupId converter
        """

        geneinfo_db = self._import_geneinfo_db(genus)
        dict_VF = dict()
        for VFgene_name, locTs in self.dict_vir.items():
            VFclass = self.vir2class[VFgene_name]

            #convert locTs to tID
            tIDs = [tID for locT in locTs for tID in geneinfo_db if locT == geneinfo_db[tID]['locT']]
            gIDs = set([self.groupmember_db.get(tID) for tID in tIDs if self.groupmember_db.get(tID)])
            for gID in gIDs:

                if gID not in dict_VF:
                    dict_VF[gID] = {'VFgene_names': VFgene_name, 'VFclass': VFclass}
                else:
                    dict_VF[gID]['VFgene_names'] += '|{}'.format(VFgene_name)
                    dict_VF[gID]['VFclass'] += '|{}'.format(VFclass)
            
        return dict_VF


    def get_ontology(self):
        """
        Check if GO ontology file exists and download if neccessary
        :return <str>: Path to GO ontology file
        """

        if os.path.isfile(self.obo_file) == False:
            print('WARNING: GO ontology not found - downloading required file')
            self._download_GO_ontology()
        return self.obo_file

    
    def _download_GO_ontology(self):
        """
        Download GO basic ontology file in obo file format
        :param obo_file <str>: Path to obo file
        """

        url = "http://purl.obolibrary.org/obo/go/go-basic.obo"
        urllib.request.urlretrieve(url, self.obo_file)
        print('    basic GO version downloaded (go-basic.obo)')


    def get_GO_terms(self, return_dict=True):
        """
        :return <dict>: db of GO terms
        """

        if not self.goterm_db:
            with gzip.GzipFile(self.goterm_db_file, 'r') as db:
                self.goterm_db = json.loads(db.read().decode('utf-8'))
        
        if return_dict == True:
            return self.goterm_db
