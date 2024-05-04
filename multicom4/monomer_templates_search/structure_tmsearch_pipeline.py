import copy
import os
import sys
import time
from multicom4.common.util import makedir_if_not_exists, check_dirs
import pandas as pd
from multicom4.monomer_templates_concatenation.parsers import TemplateHit
from multicom4.monomer_structure_refinement.iterative_refine_pipeline_v4_50 import assess_foldseek_hit, PrefilterError
from multicom4.tool import tmsearch
import dataclasses
from multicom4.monomer_structure_refinement.util import convert_taln_seq_to_a3m, build_alignment_indices
import datetime

class monomer_tmsearch_based_template_search_pipeline:

    def __init__(self, params):

        self.params = params

        self.tmsearch_database = params['tmsearch_database']

        self.template_searcher = tmsearch.TMSearch(
            binary_path=params['tmsearch_binary'],
            program_path=params['tmsearch_program_dir'],
            tmsearch_database=self.tmsearch_database)
        
        release_date_df = pd.read_csv(params['pdb_release_date_file'])
        self._release_dates = dict(zip(release_date_df['pdbcode'], release_date_df['release_date']))
        self._max_template_date = datetime.datetime.strptime(params['max_template_date'], '%Y-%m-%d')

    def create_template_df(self, templates, query_sequence):
        row_list = []
        for i in range(len(templates)):
            query = templates.loc[i, 'query']
            target = templates.loc[i, 'target'].replace('.pdb', '')
            qaln = templates.loc[i, 'qaln']
            qstart = int(templates.loc[i, 'qstart'])
            qend = int(templates.loc[i, 'qend'])
            taln = templates.loc[i, 'taln']
            tstart = templates.loc[i, 'tstart']
            tend = templates.loc[i, 'tend']
            evalue = float(templates.loc[i, 'evalue'])
            alnlen = int(templates.loc[i, 'alnlen'])

            pdbcode = target[0:4]
            if pdbcode.lower()[:4] in self._release_dates:
                hit_release_date = datetime.datetime.strptime(self._release_dates[target.lower()[:4]], '%Y-%m-%d')
                if hit_release_date < self._max_template_date:
                    hit = TemplateHit(index=i,
                                      name=target.split('.')[0],
                                      aligned_cols=alnlen,
                                      query=qaln,
                                      hit_sequence=taln,
                                      indices_query=build_alignment_indices(qaln, qstart),
                                      indices_hit=build_alignment_indices(taln, tstart),
                                      sum_probs=0.0)

                    try:
                        assess_foldseek_hit(hit=hit, query_sequence=query_sequence)
                    except PrefilterError as e:
                        msg = f'hit {hit.name.split()[0]} did not pass prefilter: {str(e)}'
                        print(msg)
                        continue

                    row_dict = dict(query=query,
                                    target=target.split()[0],
                                    taln=taln,
                                    tstart=tstart,
                                    tend=tend,
                                    qaln=qaln,
                                    qstart=qstart,
                                    qend=qend,
                                    evalue=evalue,
                                    alnlen=alnlen)

                    row_list += [row_dict]

        if len(row_list) == 0:
            return pd.DataFrame(columns=['query', 'target', 'qaln', 'taln', 'qstart', 'qend', 
                                         'tstart', 'tend', 'evalue', 'alnlen'])

        return pd.DataFrame(row_list)

    def search(self, sequence, inpdb, outdir, max_template_count = 50):
        
        # sequence = open(fasta_path).readlines()[1].rstrip('\n').lstrip('>')

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
            template = pdb_templates_result_df_filtered.loc[i, 'target']
            os.system(f"cp {self.tmsearch_database}/{template}.pdb {template_dir}")
            copied_count += 1
            if copied_count >= max_template_count:
                break

        return outfile, template_dir
