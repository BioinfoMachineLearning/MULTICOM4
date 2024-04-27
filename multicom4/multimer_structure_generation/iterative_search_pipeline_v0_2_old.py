import copy
import os
import sys
import time
from multicom4.common.util import makedir_if_not_exists, check_dirs
import pandas as pd
from multiprocessing import Pool
import dataclasses
from multicom4.tool.foldseek import *
import pickle
import numpy as np
from multicom4.monomer_templates_concatenation.parsers import TemplateHit
# from multicom4.multimer_structure_refinement.iterative_refine_pipeline_heteromer_v1_with_monomer import *
from multicom4.multimer_structure_refinement.util import convert_taln_seq_to_a3m, \
        check_and_rank_monomer_templates_local_and_global, create_template_df_with_index, \
        combine_a3ms, assess_complex_templates_homo, assess_complex_templates
from multicom4.monomer_alignment_generation.alignment import read_a3m
from multicom4.common.protein import complete_result
from multicom4.common import config

class Multimer_iterative_generation_pipeline_monomer_old(config.pipeline):

    def __init__(self, params, max_template_count=50, is_homomers=False, config_name=""):
        
        is_human = True if params['is_human'] == "1" else False

        super().__init__(is_human=is_human)
        
        self.params = params

        self.max_template_count = max_template_count

        if is_homomers:
            self.predictor_config = self.homomer_config.predictors[config_name]
        else:
            self.predictor_config = self.heteromer_config.predictors[config_name]

        release_date_df = pd.read_csv(params['pdb_release_date_file'])
        self._release_dates = dict(zip(release_date_df['pdbcode'], release_date_df['release_date']))
        self._max_template_date = datetime.datetime.strptime(params['max_template_date'], '%Y-%m-%d')
    
    def search_templates_foldseek(self, inpdb, outdir):
        makedir_if_not_exists(outdir)
        foldseek_program = self.params['foldseek_program']
        foldseek_pdb_database = ""

        other_databases = []
        if self.predictor_config.foldseek_database == "esm_atlas":
            other_databases = [os.path.join(self.params['foldseek_esm_atlas_database'], database) 
                                for database in sorted(os.listdir(self.params['foldseek_esm_atlas_database'])) 
                                if database.endswith('DB')]
            print(f"Total {len(other_databases)} to be searched!")
        else:
            foldseek_pdb_database = self.params['foldseek_pdb_database']
            alphafolddb_databases = [os.path.join(self.params['foldseek_af_database'], database) 
                                    for database in sorted(os.listdir(self.params['foldseek_af_database'])) 
                                    if database.endswith('DB')]

            other_databases += alphafolddb_databases

        foldseek_runner = Foldseek(binary_path=foldseek_program, pdb_database=foldseek_pdb_database,
                                   max_template_date=self._max_template_date, release_dates=self._release_dates,
                                   other_databases=other_databases)

        multiprocess = False
        if len(other_databases) >= 10:
            multiprocess = True

        return foldseek_runner.query(pdb=inpdb, outdir=outdir, progressive_threshold=2000, multiprocess=multiprocess)

    # def concatenate_msa_and_templates(self,
    #                                   chain_id_map,
    #                                   template_results,
    #                                   alphafold_monomer_a3ms,
    #                                   outpath):

    #     prev_df = None
    #     for i, chain_id in enumerate(chain_id_map):
    #         templates = template_results[i]['all_alignment']
    #         curr_df = create_template_df(templates)
    #         curr_df = curr_df.add_suffix(f"{i + 1}")
    #         curr_df['tpdbcode'] = curr_df[f'tpdbcode{i + 1}']
    #         curr_df = curr_df.drop([f'tpdbcode{i + 1}'], axis=1)
    #         if prev_df is None:
    #             prev_df = curr_df
    #         else:
    #             prev_df = prev_df.merge(curr_df, how="inner", on='tpdbcode')

    #     keep_indices = []
    #     chain_template_multimer_msas = {}
    #     for chain_id in chain_id_map:
    #         chain_template_multimer_msas[chain_id] = {'desc': [chain_id],
    #                                                   'seq': [chain_id_map[chain_id].sequence]}

    #     print(prev_df)
    #     for i in range(len(prev_df)):
    #         template_infos = []
    #         for j, chain_id in enumerate(chain_id_map):
    #             template = prev_df.loc[i, f'template{j + 1}']
    #             qaln = prev_df.loc[i, f'aln_query{j + 1}']
    #             qstart = int(prev_df.loc[i, f'qstart{j + 1}'])
    #             qend = int(prev_df.loc[i, f'qend{j + 1}'])
    #             taln = prev_df.loc[i, f'aln_temp{j + 1}']
    #             tstart = prev_df.loc[i, f'tstart{j + 1}']
    #             tend = prev_df.loc[i, f'tend{j + 1}']
    #             evalue = float(prev_df.loc[i, f'evalue{j + 1}'])
    #             row_dict = dict(chainid=chain_id,
    #                             template=template,
    #                             tpdbcode=template[0:4],
    #                             aln_temp=taln,
    #                             tstart=tstart,
    #                             tend=tend,
    #                             aln_query=qaln,
    #                             qstart=qstart,
    #                             qend=qend,
    #                             evalue=evalue)
    #             template_infos += [row_dict]

    #         if not assess_complex_templates(chain_id_map, template_infos):
    #             continue

    #         keep_indices += [i]
    #         for j, chain_id in enumerate(chain_id_map):
    #             query_non_gaps = [res != '-' for res in prev_df.loc[i, f'aln_query{j + 1}']]
    #             out_sequence = ''.join(convert_taln_seq_to_a3m(query_non_gaps, prev_df.loc[i, f'aln_temp{j + 1}']))
    #             aln_full = ['-'] * len(chain_id_map[chain_id].sequence)

    #             qstart = int(prev_df.loc[i, f'qstart{j + 1}'])
    #             qend = int(prev_df.loc[i, f'qend{j + 1}'])
    #             aln_full[qstart - 1:qend] = out_sequence
    #             taln_full_seq = ''.join(aln_full)
    #             chain_template_multimer_msas[chain_id]['desc'] += [prev_df.loc[i, f'template{j + 1}']]
    #             chain_template_multimer_msas[chain_id]['seq'] += [taln_full_seq]

    #     msa_out_path = outpath
    #     makedir_if_not_exists(msa_out_path)

    #     out_multimer_msas = []
    #     out_monomer_msas = []
    #     for chain_idx, chain_id in enumerate(chain_id_map):
    #         fasta_chunks = (f">{chain_template_multimer_msas[chain_id]['desc'][i]}\n" \
    #                         f"{chain_template_multimer_msas[chain_id]['seq'][i]}"
    #                         for i in range(len(chain_template_multimer_msas[chain_id]['desc'])))
    #         chain_msa_temp_interact = os.path.join(msa_out_path, chain_id + '.temp.interact')
    #         with open(chain_msa_temp_interact, 'w') as fw:
    #             fw.write('\n'.join(fasta_chunks) + '\n')

    #         out_multimer_msa = os.path.join(msa_out_path, chain_id + '.iteration.multimer.a3m')
    #         os.system(f"cp {chain_msa_temp_interact} {out_multimer_msa}")

    #         out_multimer_msas += [out_multimer_msa]

    #         monomer_template_msas = {'desc': [], 'seq': []}
    #         seen_seqs = [chain_template_multimer_msas[chain_id]['seq'][i]
    #                      for i in range(len(chain_template_multimer_msas[chain_id]['desc']))]
    #         templates = template_results[chain_idx]['all_alignment']
    #         for i in range(len(templates)):
    #             query_non_gaps = [res != '-' for res in templates.loc[i, f'qaln']]
    #             out_sequence = ''.join(convert_taln_seq_to_a3m(query_non_gaps, templates.loc[i, f'taln']))
    #             aln_full = ['-'] * len(chain_id_map[chain_id].sequence)

    #             qstart = int(templates.loc[i, f'qstart'])
    #             qend = int(templates.loc[i, f'qend'])
    #             aln_full[qstart - 1:qend] = out_sequence
    #             taln_full_seq = ''.join(aln_full)

    #             if taln_full_seq not in seen_seqs:
    #                 monomer_template_msas['desc'] += [templates.loc[i, 'target']]
    #                 monomer_template_msas['seq'] += [taln_full_seq]

    #         fasta_chunks = (f">{monomer_template_msas['desc'][i]}\n{monomer_template_msas['seq'][i]}"
    #                         for i in range(len(monomer_template_msas['desc'])))
            
    #         msa_temp_monomer = os.path.join(msa_out_path, chain_id + '.temp.monomer')
    #         with open(msa_temp_monomer, 'w') as fw:
    #             fw.write('\n'.join(fasta_chunks) + '\n')

    #         iteration_monomer_a3m = os.path.join(msa_out_path, chain_id + ".iteration.monomer.a3m")
    #         combine_a3ms([alphafold_monomer_a3ms[chain_idx], msa_temp_monomer], iteration_monomer_a3m)
    #         out_monomer_msas += [iteration_monomer_a3m]

    #     interact_dict = {}
    #     msa_len = -1
    #     for i in range(0, len(out_multimer_msas)):
    #         msa_sequences, msa_descriptions = parse_fasta(out_multimer_msas[i])
    #         current_len = len(msa_descriptions)
    #         if msa_len == -1:
    #             msa_len = current_len
    #         elif current_len != msa_len:
    #             raise Exception(f"The length of each msas are not equal! {out_multimer_msas}")
    #         interact_dict[f'index_{i + 1}'] = [j for j in range(msa_len)]

    #     interact_df = pd.DataFrame(interact_dict)
    #     interact_csv = os.path.join(outpath, 'interaction.iteration.csv')
    #     interact_df.to_csv(interact_csv)

    #     top_template_files = []
    #     for template_result, chain_id in zip(template_results, chain_id_map):
    #         top_template_file = os.path.join(outpath, f"{chain_id}.top{self.max_template_count}")
    #         check_and_rank_monomer_templates_local_and_global(template_result=template_result,
    #                                                           outfile=top_template_file,
    #                                                           query_sequence=chain_id_map[chain_id].sequence,
    #                                                           max_template_count=self.max_template_count,
    #                                                           release_dates=self._release_dates,
    #                                                           max_template_date=self._max_template_date)
    #         top_template_files += [top_template_file]

    #     return top_template_files, out_multimer_msas, out_monomer_msas, interact_csv

    def copy_atoms_and_unzip(self, templates, outdir):
        os.makedirs(outdir, exist_ok=True)
        os.chdir(outdir)
        num_templates = min(len(templates), 50)
        for i in range(num_templates):
            template_pdb = templates.loc[i, 'target']
            trg_pdb_path = os.path.join(outdir, template_pdb)
            if os.path.exists(trg_pdb_path):
                continue
            if template_pdb.find('.pdb') > 0:
                template_path = find_template_in_alphafolddb(self.params['foldseek_af_database_dir'], template_pdb)
                if template_path is None:
                    http_address = "https://sysbio.rnet.missouri.edu/multicom_cluster/multicom3_db_tools/databases/af_pdbs"
                    os.system(f"wget -P {outdir} --no-check-certificate {http_address}/{template_pdb}")
                else:
                    os.system(f"cp {template_path} {outdir}")

            else:
                template_path = os.path.join(self.params['foldseek_pdb_database_dir'], template_pdb)
                os.system(f"cp {template_path} {outdir}")
                os.system(f"gunzip -f {template_pdb}")

    # def search_single(self, fasta_file, chain_id_map, monomer_pdb_dirs, monomer_alphafold_a3ms, outdir):

    #     fasta_file = os.path.abspath(fasta_file)

    #     targetname = pathlib.Path(fasta_file).stem

    #     print(f"Processing {targetname}")

    #     outdir = os.path.abspath(outdir)

    #     makedir_if_not_exists(outdir)

    #     cwd = os.getcwd()

    #     makedir_if_not_exists(outdir)

    #     out_model_dir = outdir

    #     prepare_dir = os.path.join(outdir, 'prepare')

    #     makedir_if_not_exists(prepare_dir)

    #     common_parameters =   f"--fasta_path={fasta_file} " \
    #                           f"--env_dir={self.params['alphafold_env_dir']} " \
    #                           f"--database_dir={self.params['alphafold_database_dir']} " \
    #                           f"--benchmark={self.params['alphafold_benchmark']} " \
    #                           f"--use_gpu_relax={self.params['use_gpu_relax']} " \
    #                           f"--max_template_date={self.params['max_template_date']} "

    #     multimer_num_ensemble = self.get_heteromer_config(self.predictor_config, 'num_ensemble')
    #     multimer_num_recycle = self.get_heteromer_config(self.predictor_config, 'num_recycle')
    #     num_multimer_predictions_per_model = self.get_heteromer_config(self.predictor_config, 'predictions_per_model')
    #     model_preset = self.get_heteromer_config(self.predictor_config, 'model_preset')
    #     relax_topn_predictions = self.get_heteromer_config(self.predictor_config, 'relax_topn_predictions')
    #     dropout = self.get_heteromer_config(self.predictor_config, 'dropout')
    #     dropout_structure_module = self.get_heteromer_config(self.predictor_config, 'dropout_structure_module')     

    #     common_parameters +=  f"--multimer_num_ensemble={multimer_num_ensemble} " \
    #                             f"--multimer_num_recycle={multimer_num_recycle} " \
    #                             f"--num_multimer_predictions_per_model={num_multimer_predictions_per_model} " \
    #                             f"--model_preset={model_preset} " \
    #                             f"--relax_topn_predictions={relax_topn_predictions} " \
    #                             f"--models_to_relax=TOPN "


    #     if not complete_result(out_model_dir, 5 * num_multimer_predictions_per_model):

    #         out_template_dir = os.path.join(outdir, 'templates')

    #         makedir_if_not_exists(out_template_dir)

    #         template_results = []
    #         alphafold_monomer_a3ms = []

    #         for chain_id in chain_id_map:

    #             monomer_work_dir = os.path.join(prepare_dir, chain_id)

    #             makedir_if_not_exists(monomer_work_dir)

    #             if not os.path.exists(monomer_alphafold_a3ms[chain_id]):
    #                 raise Exception(f"Cannot find the monomer final a3m in {monomer_alphafold_a3ms[chain_id]}")

    #             alphafold_monomer_a3m = os.path.join(outdir, chain_id + ".alphafold.monomer.a3m")
    #             os.system(f"cp {monomer_alphafold_a3ms[chain_id]} {alphafold_monomer_a3m}")

    #             alphafold_monomer_a3ms += [alphafold_monomer_a3m]

    #             chain_pdb = os.path.join(monomer_work_dir, chain_id + '.pdb')

    #             os.system(f"cp {monomer_pdb_dirs[chain_id]} {chain_pdb}")

    #             foldseek_res = self.search_templates_foldseek(inpdb=chain_pdb, outdir=os.path.join(monomer_work_dir, 'foldseek'))

    #             if len(foldseek_res['all_alignment']) == 0:
    #                 print(
    #                     f"Cannot find any templates for {chain_id}")
    #                 break

    #             template_results += [foldseek_res]

    #             self.copy_atoms_and_unzip(templates=foldseek_res['all_alignment'],
    #                                       outdir=out_template_dir)

    #         if len(template_results) != len(chain_id_map):
    #             return None

    #         template_files, multimer_msa_files, monomer_msa_files, msa_pair_file = \
    #             self.concatenate_msa_and_templates(chain_id_map=chain_id_map,
    #                                                template_results=template_results,
    #                                                alphafold_monomer_a3ms=alphafold_monomer_a3ms,
    #                                                outpath=outdir)

    #         cmd = f"python {self.params['alphafold_multimer_program']} " + common_parameters + f"--output_dir={out_model_dir} "

    #         cmd += f"--monomer_a3ms={','.join(monomer_msa_files)} "
    #         if self.predictor_config.msa_paired_source == "None":
    #             cmd += "--msa_pairing_hetero=false "
    #         else:
    #             cmd += f"--multimer_a3ms={','.join(multimer_msa_files)} " \
    #                    f"--msa_pair_file={msa_pair_file} " \

    #         if self.predictor_config.template_source == "notemplate":
    #             cmd += "--notemplate=true "
    #         else:
    #             cmd +=  f"--monomer_temp_csvs={','.join(template_files)} " \
    #                     f"--struct_atom_dir={out_template_dir} "

    #         try:
    #             os.chdir(self.params['alphafold_program_dir'])
    #             print(cmd)
    #             os.system(cmd)
    #         except Exception as e:
    #             print(e)

    #     os.chdir(cwd)

    def concatenate_msa_and_templates_homo(self,
                                           chain_id_map,
                                           template_results,
                                           alphafold_monomer_a3ms,
                                           outpath):

        prev_df = None
        for i, chain_id in enumerate(chain_id_map):
            templates = template_results[i]['all_alignment']
            curr_df = create_template_df_with_index(templates)
            curr_df = curr_df.add_suffix(f"{i + 1}")
            curr_df['index'] = curr_df[f'index{i + 1}']
            curr_df = curr_df.drop([f'index{i + 1}'], axis=1)
            if prev_df is None:
                prev_df = curr_df
            else:
                prev_df = prev_df.merge(curr_df, how="inner", on='index')

        keep_indices = []
        chain_template_multimer_msas = {}
        for chain_id in chain_id_map:
            chain_template_multimer_msas[chain_id] = {'desc': [chain_id],
                                                      'seq': [chain_id_map[chain_id].sequence]}

        print(prev_df)
        for i in range(len(prev_df)):
            template_infos = []
            for j, chain_id in enumerate(chain_id_map):
                template = prev_df.loc[i, f'template{j + 1}']
                qaln = prev_df.loc[i, f'aln_query{j + 1}']
                qstart = int(prev_df.loc[i, f'qstart{j + 1}'])
                qend = int(prev_df.loc[i, f'qend{j + 1}'])
                taln = prev_df.loc[i, f'aln_temp{j + 1}']
                tstart = prev_df.loc[i, f'tstart{j + 1}']
                tend = prev_df.loc[i, f'tend{j + 1}']
                evalue = float(prev_df.loc[i, f'evalue{j + 1}'])
                row_dict = dict(chainid=chain_id,
                                template=template,
                                tpdbcode=template[0:4],
                                aln_temp=taln,
                                tstart=tstart,
                                tend=tend,
                                aln_query=qaln,
                                qstart=qstart,
                                qend=qend,
                                evalue=evalue)
                template_infos += [row_dict]

            # if not assess_complex_templates(chain_id_map, template_infos):
            #     continue

            keep_indices += [i]
            for j, chain_id in enumerate(chain_id_map):
                query_non_gaps = [res != '-' for res in prev_df.loc[i, f'aln_query{j + 1}']]
                out_sequence = ''.join(convert_taln_seq_to_a3m(query_non_gaps, prev_df.loc[i, f'aln_temp{j + 1}']))
                aln_full = ['-'] * len(chain_id_map[chain_id].sequence)

                qstart = int(prev_df.loc[i, f'qstart{j + 1}'])
                qend = int(prev_df.loc[i, f'qend{j + 1}'])
                aln_full[qstart - 1:qend] = out_sequence
                taln_full_seq = ''.join(aln_full)
                chain_template_multimer_msas[chain_id]['desc'] += [prev_df.loc[i, f'template{j + 1}']]
                chain_template_multimer_msas[chain_id]['seq'] += [taln_full_seq]

        msa_out_path = outpath
        makedir_if_not_exists(msa_out_path)

        out_monomer_msas = []
        for chain_idx, chain_id in enumerate(chain_id_map):
            fasta_chunks = (f">{chain_template_multimer_msas[chain_id]['desc'][i]}\n" \
                            f"{chain_template_multimer_msas[chain_id]['seq'][i]}"
                            for i in range(len(chain_template_multimer_msas[chain_id]['desc'])))
            
            msa_temp_interact = os.path.join(msa_out_path, chain_id + '.temp.interact')
            with open(msa_temp_interact, 'w') as fw:
                fw.write('\n'.join(fasta_chunks) + '\n')

            iteration_multimer_a3m = os.path.join(msa_out_path, chain_id + '.iteration.multimer.a3m')
            os.system(f"cp {msa_temp_interact} {iteration_multimer_a3m}")

            monomer_template_msas = {'desc': [], 'seq': []}
            seen_seqs = [chain_template_multimer_msas[chain_id]['seq'][i]
                         for i in range(len(chain_template_multimer_msas[chain_id]['desc']))]
            templates = template_results[chain_idx]['all_alignment']
            for i in range(len(templates)):
                query_non_gaps = [res != '-' for res in templates.loc[i, f'qaln']]
                out_sequence = ''.join(convert_taln_seq_to_a3m(query_non_gaps, templates.loc[i, f'taln']))
                aln_full = ['-'] * len(chain_id_map[chain_id].sequence)

                qstart = int(templates.loc[i, f'qstart'])
                qend = int(templates.loc[i, f'qend'])
                aln_full[qstart - 1:qend] = out_sequence
                taln_full_seq = ''.join(aln_full)

                if taln_full_seq not in seen_seqs:
                    monomer_template_msas['desc'] += [templates.loc[i, 'target']]
                    monomer_template_msas['seq'] += [taln_full_seq]

            fasta_chunks = (f">{monomer_template_msas['desc'][i]}\n{monomer_template_msas['seq'][i]}"
                            for i in range(len(monomer_template_msas['desc'])))
            msa_temp_monomer = os.path.join(msa_out_path, chain_id + '.temp.monomer')
            with open(msa_temp_monomer, 'w') as fw:
                fw.write('\n'.join(fasta_chunks) + '\n')

            iteration_monomer_a3m = os.path.join(msa_out_path, chain_id + '.iteration.monomer.a3m')
            combine_a3ms([iteration_multimer_a3m,
                          alphafold_monomer_a3ms[chain_idx],
                          msa_temp_monomer],
                        iteration_monomer_a3m)
            out_monomer_msas += [iteration_monomer_a3m]

        top_template_files = []
        for template_result, chain_id in zip(template_results, chain_id_map):
            top_template_file = os.path.join(outpath, f"{chain_id}.top{self.max_template_count}")
            check_and_rank_monomer_templates_local_and_global(template_result=template_result,
                                                              outfile=top_template_file,
                                                              query_sequence=chain_id_map[chain_id].sequence,
                                                              max_template_count=self.max_template_count)
            top_template_files += [top_template_file]

        return top_template_files, out_monomer_msas

    def search_single_homo(self, fasta_file, chain_id_map, monomer_pdb_dirs, 
                           monomer_alphafold_a3ms, outdir, monomer_template_stos = []):

        fasta_file = os.path.abspath(fasta_file)

        targetname = pathlib.Path(fasta_file).stem

        print(f"Processing {targetname}")

        outdir = os.path.abspath(outdir)

        makedir_if_not_exists(outdir)

        cwd = os.getcwd()

        makedir_if_not_exists(outdir)

        out_model_dir = outdir

        prepare_dir = os.path.join(outdir, 'prepare')

        makedir_if_not_exists(prepare_dir)

        common_parameters =   f"--fasta_path={fasta_file} " \
                              f"--env_dir={self.params['alphafold_env_dir']} " \
                              f"--database_dir={self.params['alphafold_database_dir']} " \
                              f"--benchmark={self.params['alphafold_benchmark']} " \
                              f"--use_gpu_relax={self.params['use_gpu_relax']} " \
                              f"--max_template_date={self.params['max_template_date']} "

        multimer_num_ensemble = self.get_homomer_config(self.predictor_config, 'num_ensemble')
        multimer_num_recycle = self.get_homomer_config(self.predictor_config, 'num_recycle')
        num_multimer_predictions_per_model = self.get_homomer_config(self.predictor_config, 'predictions_per_model')
        model_preset = self.get_homomer_config(self.predictor_config, 'model_preset')
        relax_topn_predictions = self.get_homomer_config(self.predictor_config, 'relax_topn_predictions')
        dropout = self.get_homomer_config(self.predictor_config, 'dropout')
        dropout_structure_module = self.get_homomer_config(self.predictor_config, 'dropout_structure_module')     

        common_parameters +=  f"--multimer_num_ensemble={multimer_num_ensemble} " \
                                f"--multimer_num_recycle={multimer_num_recycle} " \
                                f"--num_multimer_predictions_per_model={num_multimer_predictions_per_model} " \
                                f"--model_preset={model_preset} " \
                                f"--relax_topn_predictions={relax_topn_predictions} " \
                                f"--models_to_relax=TOPN "

        if not complete_result(out_model_dir, 5 * num_multimer_predictions_per_model):

            out_template_dir = os.path.join(prepare_dir, 'templates')

            template_results = []
            alphafold_monomer_a3ms = []

            for chain_id in chain_id_map:

                monomer_work_dir = os.path.join(prepare_dir, chain_id)

                makedir_if_not_exists(monomer_work_dir)

                if not os.path.exists(monomer_alphafold_a3ms[chain_id]):
                    raise Exception(f"Cannot find the monomer final a3m in {monomer_alphafold_a3ms[chain_id]}")

                monomer_alphafold_a3m = os.path.join(outdir, chain_id + ".alphafold.monomer.a3m")
                os.system(f"cp {monomer_alphafold_a3ms[chain_id]} {monomer_alphafold_a3m}")

                alphafold_monomer_a3ms += [monomer_alphafold_a3m]

                chain_pdb = os.path.join(monomer_work_dir, chain_id + '.pdb')

                os.system(f"cp {monomer_pdb_dirs[chain_id]} {chain_pdb}")

                foldseek_res = self.search_templates_foldseek(inpdb=chain_pdb, outdir=os.path.join(monomer_work_dir, 'foldseek'))

                if len(foldseek_res['all_alignment']) == 0:
                    print(f"Cannot find any templates for {chain_id}")
                    break

                template_results += [foldseek_res]

            if len(template_results) != len(chain_id_map):
                return None

            template_files, monomer_msa_files = \
                self.concatenate_msa_and_templates_homo(chain_id_map=chain_id_map,
                                                        template_results=template_results,
                                                        alphafold_monomer_a3ms=alphafold_monomer_a3ms,
                                                        outpath=outdir)

            if self.predictor_config.template_source == "foldseek":
                for template_file in template_files:
                    templates = pd.read_csv(template_file, sep='\t')
                    self.copy_atoms_and_unzip(templates=templates, outdir=out_template_dir)

            cmd = f"python {self.params['alphafold_multimer_program']} " + common_parameters + f"--output_dir={out_model_dir} "

            cmd += f"--monomer_a3ms={','.join(monomer_msa_files)} "

            if self.predictor_config.template_source == "notemplate":
                cmd += "--notemplate=true "
            elif self.predictor_config.template_source == "pdb70":
                cmd += f"--template_stos {','.join(template_stos)} "
            elif self.predictor_config.template_source == "foldseek":
                cmd +=  f"--monomer_temp_csvs={','.join(template_files)} " \
                        f"--struct_atom_dir={out_template_dir} "

            try:
                os.chdir(self.params['alphafold_program_dir'])
                print(cmd)
                # os.system(cmd)
            except Exception as e:
                print(e)

        os.chdir(cwd)
