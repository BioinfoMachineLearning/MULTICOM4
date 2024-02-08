import copy
import os
import sys
import time
from multicom4.common.util import makedir_if_not_exists, check_dirs
import pandas as pd
from multiprocessing import Pool
from multicom4.monomer_templates_concatenation import parsers
from multicom4.monomer_structure_refinement.iterative_refine_pipeline_v4_50 import assess_foldseek_hit
from multicom4.tool import tmsearch
import dataclasses
from multicom4.monomer_structure_refinement.util import *

def create_df(hits):
    row_list = []
    for index, hit in enumerate(hits):
        row_dict = dict(index=index,
                        name=hit.name,
                        tpdbcode=hit.name[0:4],
                        aligned_cols=hit.aligned_cols,
                        sum_probs=hit.sum_probs,
                        query=hit.query,
                        hit_sequence=hit.hit_sequence,
                        indices_query='_'.join([str(i) for i in hit.indices_query]),
                        indices_hit='_'.join([str(i) for i in hit.indices_hit]))
        row_list += [row_dict]

    if len(row_list) == 0:
        empty_dict = dict(index=0,
                          name='empty',
                          tpdbcode='empty',
                          aligned_cols=0,
                          sum_probs=0,
                          query='empty',
                          hit_sequence='empty',
                          indices_query=None,
                          indices_hit=None)
        row_list += [empty_dict]
    return pd.DataFrame(row_list)


class monomer_tmsearch_based_template_search_pipeline:

    def __init__(self, params):

        self.params = params

        self.tmsearch_database = params['tmsearch_database']

        self.template_searcher = tmsearch.TMSearch(
            binary_path=params['tmsearch_binary'],
            program_path=params['tmsearch_program_dir'],
            tmsearch_database=self.tmsearch_database)
        
        # release_date_df = pd.read_csv(params['pdb_release_date_file'])
        # self._release_dates = dict(zip(release_date_df['pdbcode'], release_date_df['release_date']))
        # self._max_template_date = datetime.datetime.strptime(params['max_template_date'], '%Y-%m-%d')

    def create_template_df(self, templates, query_sequence, check_hit=True):
        row_list = []
        for i in range(len(templates)):
            target = templates.loc[i, 'target']
            qaln = templates.loc[i, 'qaln']
            qstart = int(templates.loc[i, 'qstart'])
            qend = int(templates.loc[i, 'qend'])
            taln = templates.loc[i, 'taln']
            tstart = templates.loc[i, 'tstart']
            tend = templates.loc[i, 'tend']
            evalue = float(templates.loc[i, 'evalue'])
            aln_len = int(templates.loc[i, 'alnlen'])

            pdbcode = target[0:4]

            hit = parsers.TemplateHit(index=i,
                                    name=templates.loc[i, 'target'].split('.')[0],
                                    aligned_cols=int(templates.loc[i, 'alnlen']),
                                    query=templates.loc[i, 'qaln'],
                                    hit_sequence=templates.loc[i, 'taln'],
                                    indices_query=build_alignment_indices(templates.loc[i, 'qaln'],
                                                                            templates.loc[i, 'qstart']),
                                    indices_hit=build_alignment_indices(templates.loc[i, 'taln'],
                                                                        templates.loc[i, 'tstart']),
                                    sum_probs=0.0)

            if check_hit:
                try:
                    assess_foldseek_hit(hit=hit, query_sequence=query_sequence)
                except PrefilterError as e:
                    msg = f'hit {hit.name.split()[0]} did not pass prefilter: {str(e)}'
                    print(msg)
                    continue

            row_dict = dict(template=target.split()[0],
                            tpdbcode=pdbcode,
                            aln_temp=taln,
                            tstart=tstart,
                            tend=tend,
                            aln_query=qaln,
                            qstart=qstart,
                            qend=qend,
                            tmscore=evalue,
                            aligned_length=aln_len)

            row_list += [row_dict]
        if len(row_list) == 0:
            return pd.DataFrame(columns=['template', 'tpdbcode', 'aln_temp', 'tstart', 'tend',
                                        'aln_query', 'qstart', 'qend', 'evalue', 'aligned_length'])
        return pd.DataFrame(row_list)

    def search(self, fasta_path, inpdb, outdir, max_template_count = 50):
        
        sequence = open(fasta_path).readlines()[1].rstrip('\n').lstrip('>')

        makedir_if_not_exists(outdir)
        outfile = os.path.join(outdir, 'tmsearch_templates.csv')
        template_dir = os.path.join(outdir, 'templates')
        makedir_if_not_exists(template_dir)

        pdb_templates_result_df = self.template_searcher.query(inpdb, outdir)

        pdb_templates_result_df_filtered = self.create_template_df(pdb_templates_result_df, sequence)
        print(pdb_templates_result_df_filtered)
        pdb_templates_result_df_filtered = pdb_templates_result_df_filtered.sort_values(by='evalue', ascending=False)
        pdb_templates_result_df_filtered.reset_index(inplace=True, drop=True)
        pdb_templates_result_df_filtered.to_csv(outfile, sep='\t') 

        copied_count = 0
        for i in range(len(pdb_templates_result_df_filtered)):
            template = pdb_templates_result_df_filtered.iloc[i, 'template']
            os.system(f"cp {self.tmsearch_database}/{template}.pdb {template_dir}")

        return outfile, template_dir
