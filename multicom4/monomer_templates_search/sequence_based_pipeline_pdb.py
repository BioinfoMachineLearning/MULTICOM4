import copy
import os
import sys
import time
from multicom4.common.util import makedir_if_not_exists, check_dirs
import pandas as pd
from multiprocessing import Pool
from multicom4.monomer_templates_concatenation import parsers
from multicom4.tool import hhsearch
from multicom4.tool import hhalign
import dataclasses
import datetime
from typing import Any, Dict, Mapping, Optional, Sequence, Tuple

# Prefilter exceptions.
class PrefilterError(Exception):
    """A base class for template prefilter exceptions."""


class DateError(PrefilterError):
    """An error indicating that the hit date was after the max allowed date."""


class AlignRatioError(PrefilterError):
    """An error indicating that the hit align ratio to the query was too small."""


class DuplicateError(PrefilterError):
    """An error indicating that the hit was an exact subsequence of the query."""


class LengthError(PrefilterError):
    """An error indicating that the hit was too short."""


def create_df(targetname, hits):
    row_list = []
    for index, hit in enumerate(hits):
        row_dict = dict(query=targetname,
                        target=hit.name,
                        alnlen=hit.aligned_cols,
                        sum_probs=hit.sum_probs,
                        qaln=hit.query,
                        qstart=hit.indices_query[0]+1,
                        qend=hit.indices_query[len(hit.indices_query)-1]+1,
                        taln=hit.hit_sequence,
                        tstart=hit.indices_hit[0] + 1,
                        tend=hit.indices_hit[len(hit.indices_hit) - 1] + 1)
        row_list += [row_dict]

    if len(row_list) == 0:
        return pd.DataFrame(columns=['query', 'target', 'alnlen', 'sum_probs', 'qaln', 'qstart', 'qend',
                                     'taln', 'tstart', 'tend'])
    return pd.DataFrame(row_list)


def assess_hhsearch_hit(
        hit: parsers.TemplateHit,
        query_sequence: str,
        max_template_date: datetime.datetime,
        release_dates: Mapping[str, datetime.datetime],
        max_subsequence_ratio: float = 0.95,
        min_align_ratio: float = 0.1) -> bool:

    aligned_cols = hit.aligned_cols
    align_ratio = aligned_cols / len(query_sequence)

    template_sequence = hit.hit_sequence.replace('-', '')
    length_ratio = float(len(template_sequence)) / len(query_sequence)

    duplicate = (template_sequence in query_sequence and
                 length_ratio > max_subsequence_ratio)

    # if max_template_date is not None and release_dates is not None:
    #     if hit.name.lower()[:4] in release_dates:
    #         hit_release_date = datetime.datetime.strptime(release_dates[hit.name.lower()[:4]], '%Y-%m-%d')
    #         if hit_release_date > max_template_date:
    #             raise DateError(f'Date ({release_dates[hit.name.lower()[:4]]}) > max template date '
    #                             f'({max_template_date}).')
    #     else:
    #         raise DateError(f'Cannot find release date for ({hit.name.lower()[:4]}).')

    if align_ratio <= min_align_ratio:
        raise AlignRatioError('Proportion of residues aligned to query too small. '
                              f'Align ratio: {align_ratio}.')

    if duplicate:
        raise DuplicateError('Template is an exact subsequence of query with large '
                             f'coverage. Length ratio: {length_ratio}.')

    if len(template_sequence) < 10:
        raise LengthError(f'Template too short. Length: {len(template_sequence)}.')

    return True


class monomer_sequence_based_template_search_pipeline:

    def __init__(self, params):

        self.params = params

        self.template_searcher = hhsearch.HHSearch(
            binary_path=params['hhsearch_program'],
            databases=[params['pdb_sort90_hhsuite_database']],
            input_format='hmm')

        self.pdbdir = params['pdb_sort90_atom_dir']

        self.hhmake_program = params['hhmake_program']

        release_date_df = pd.read_csv(params['pdb_release_date_file'])
        self._release_dates = dict(zip(release_date_df['pdbcode'], release_date_df['release_date']))
        self._max_template_date = datetime.datetime.strptime(params['max_template_date'], '%Y-%m-%d')

    def copy_atoms_and_unzip(self, templates, outdir):
        os.chdir(outdir)
        for i in range(len(templates)):
            template_pdb = templates.loc[i, 'target']
            template_path = os.path.join(self.pdbdir, template_pdb + '.atom.gz')
            os.system(f"cp {template_path} .")
            os.system(f"gunzip -f {template_pdb}.atom.gz")

    def search(self, targetname, sequence, a3m, outdir):
        makedir_if_not_exists(outdir)
        # pdb_hits_out_path = os.path.join(outdir, f'pdb_hits.{self.template_searcher.output_format}')
        pdb_hits_out_path = os.path.join(outdir, f'output.hhr')
        if os.path.exists(pdb_hits_out_path):
            pdb_templates_result = open(pdb_hits_out_path, encoding='ISO-8859-1').read()
        else:
            trg_a3m = os.path.join(outdir, targetname + '.a3m')
            trg_hmm = os.path.join(outdir, targetname + '.hmm')
            with open(a3m) as f:
                msa_for_templates = f.read()
                if a3m.find('.sto') > 0:
                    msa_for_templates = parsers.deduplicate_stockholm_msa(msa_for_templates)
                    msa_for_templates = parsers.remove_empty_columns_from_stockholm_msa(msa_for_templates)
                    msa_for_templates = parsers.convert_stockholm_to_a3m(msa_for_templates)

                    with open(trg_a3m, 'w') as fw:
                        fw.write(msa_for_templates)
                else:
                    os.system(f"cp {a3m} {trg_a3m}")

            os.system(f"{self.hhmake_program} -i {trg_a3m} -o {trg_hmm}")
            with open(trg_hmm) as f:
                msa_for_templates = f.read()
                pdb_templates_result = self.template_searcher.query(msa_for_templates, outdir)

        pdb_template_hits = parsers.parse_hhr(hhr_string=pdb_templates_result)

        pdb_template_hits = sorted(pdb_template_hits, key=lambda x: x.sum_probs, reverse=True)

        curr_template_hits = []
        for hit in pdb_template_hits:
            try:
                assess_hhsearch_hit(hit=hit, query_sequence=sequence, max_template_date=self._max_template_date, release_dates=self._release_dates)
            except PrefilterError as e:
                msg = f'hit {hit.name.split()[0]} did not pass prefilter: {str(e)}'
                print(msg)
                continue
            curr_template_hits += [hit]

        curr_pd = create_df(targetname, curr_template_hits)

        curr_pd.to_csv(os.path.join(outdir, 'sequence_templates.csv'), sep='\t')

        template_dir = os.path.join(outdir, 'templates')

        makedir_if_not_exists(template_dir)

        self.copy_atoms_and_unzip(templates=curr_pd, outdir=template_dir)

        return os.path.join(outdir, 'sequence_templates.csv')
