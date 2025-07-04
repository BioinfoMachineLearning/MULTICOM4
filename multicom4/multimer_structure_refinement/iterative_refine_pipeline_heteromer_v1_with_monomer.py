import copy
import os
import sys
import time, json
from multicom4.common.util import makedir_if_not_exists, check_dirs
import pandas as pd
from multiprocessing import Pool
import dataclasses
from multicom4.tool.foldseek import *
import pickle
import numpy as np
from multicom4.multimer_structure_refinement.util import *
from multicom4.common.protein import complete_result
from multicom4.common import config

class Multimer_iterative_refinement_pipeline(config.pipeline):

    def __init__(self, params, config_name):
        super().__init__()
        self.params = params
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
            esm_atlas_databases = [os.path.join(self.params['foldseek_esm_atlas_database'], database) 
                                for database in sorted(os.listdir(self.params['foldseek_esm_atlas_database'])) 
                                if database.endswith('DB')]
            # e.g, tm_.80_.90_plddt_.80_.90_14.DB
            pattern = r'tm_([\d.]+)_[\d.]+_plddt_([\d.]+)_[\d.]+'
            for esm_atlas_database in esm_atlas_databases:
                match = re.search(pattern, esm_atlas_database)
                ptm = float(match.group(1))
                plddt = float(match.group(2))
                if ptm >= self.predictor_config.ptm_threshold and plddt >= self.predictor_config.plddt_threshold:
                    other_databases += [esm_atlas_database]

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

    def concatenate_msa_and_templates(self,
                                      chain_id_map,
                                      template_results,
                                      start_msa_path,
                                      outpath,
                                      iteration):

        chain_template_msas = {}
        evalue_thresholds = [1e-7, 1e-6, 1e-5, 1e-4, 1e-3]
        tmscore_thresholds = [0.8, 0.7, 0.6, 0.5, 0.4]
        complex_templates_df_filtered = None
        for evalue_threshold, tmscore_threshold in zip(evalue_thresholds, tmscore_thresholds):
            chain_template_msas = {}
            for chain_id in chain_id_map:
                chain_template_msas[chain_id] = {'desc': [chain_id],
                                                 'seq': [chain_id_map[chain_id].sequence]}

            complex_templates_df = None
            for chain_idx, (chain_id, template_result) in enumerate(zip(chain_id_map, template_results)):
                evalue_keep_indices = []
                for i in range(len(template_result['local_alignment'])):
                    if template_result['local_alignment'].loc[i, 'evalue'] < evalue_threshold:
                        evalue_keep_indices += [i]

                tmscore_keep_indices = []
                for i in range(len(template_result['global_alignment'])):
                    if template_result['global_alignment'].loc[i, 'evalue'] > tmscore_threshold:
                        tmscore_keep_indices += [i]

                templates_filtered = copy.deepcopy(template_result['local_alignment'].iloc[evalue_keep_indices])
                templates_filtered = templates_filtered.append(
                    copy.deepcopy(template_result['global_alignment'].iloc[tmscore_keep_indices]))
                templates_filtered.drop(templates_filtered.filter(regex="Unnamed"), axis=1, inplace=True)
                templates_filtered.reset_index(inplace=True, drop=True)

                curr_df = create_template_df(templates_filtered)
                # print(curr_df)
                curr_df = curr_df.add_suffix(f"{chain_idx + 1}")
                curr_df['tpdbcode'] = curr_df[f'tpdbcode{chain_idx + 1}']
                curr_df = curr_df.drop([f'tpdbcode{chain_idx + 1}'], axis=1)
                if complex_templates_df is None:
                    complex_templates_df = curr_df
                else:
                    complex_templates_df = complex_templates_df.merge(curr_df, how="inner", on='tpdbcode')

            keep_indices = []
            seen_complex_seq = []
            seen_complex_seq += ["".join([chain_template_msas[chain_id]['seq'][0] for chain_id in chain_template_msas])]
            seen_pdbcodes = []
            for i in range(len(complex_templates_df)):
                if len(keep_indices) > self.predictor_config.max_template_count:
                    break
                template_infos = []
                for j, chain_id in enumerate(chain_id_map):
                    template = complex_templates_df.loc[i, f'template{j + 1}']
                    qaln = complex_templates_df.loc[i, f'aln_query{j + 1}']
                    qstart = int(complex_templates_df.loc[i, f'qstart{j + 1}'])
                    qend = int(complex_templates_df.loc[i, f'qend{j + 1}'])
                    taln = complex_templates_df.loc[i, f'aln_temp{j + 1}']
                    tstart = complex_templates_df.loc[i, f'tstart{j + 1}']
                    tend = complex_templates_df.loc[i, f'tend{j + 1}']
                    evalue = float(complex_templates_df.loc[i, f'evalue{j + 1}'])
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

                if complex_templates_df.loc[i, 'tpdbcode'] in seen_pdbcodes:
                    continue

                if not assess_complex_templates(chain_id_map, template_infos):
                    continue

                monomer_template_seqs = {}
                for j, chain_id in enumerate(chain_id_map):
                    query_non_gaps = [res != '-' for res in complex_templates_df.loc[i, f'aln_query{j + 1}']]
                    out_sequence = ''.join(
                        convert_taln_seq_to_a3m(query_non_gaps, complex_templates_df.loc[i, f'aln_temp{j + 1}']))
                    aln_full = ['-'] * len(chain_id_map[chain_id].sequence)
                    qstart = int(complex_templates_df.loc[i, f'qstart{j + 1}'])
                    qend = int(complex_templates_df.loc[i, f'qend{j + 1}'])
                    aln_full[qstart - 1:qend] = out_sequence
                    taln_full_seq = ''.join(aln_full)
                    monomer_template_dict = {'desc': complex_templates_df.loc[i, f'template{j + 1}'],
                                             'seq': taln_full_seq}
                    monomer_template_seqs[chain_id] = monomer_template_dict

                complex_template_seq = "".join(
                    [monomer_template_seqs[chain_id]['seq'] for chain_id in monomer_template_seqs])
                if complex_template_seq not in seen_complex_seq:
                    for chainid in monomer_template_seqs:
                        chain_template_msas[chainid]['desc'] += [monomer_template_seqs[chainid]['desc']]
                        chain_template_msas[chainid]['seq'] += [monomer_template_seqs[chainid]['seq']]
                    seen_complex_seq += [complex_template_seq]
                    keep_indices += [i]
                    seen_pdbcodes += [complex_templates_df.loc[i, 'tpdbcode']]

            complex_templates_df_filtered = copy.deepcopy(complex_templates_df.iloc[keep_indices])
            complex_templates_df_filtered.drop(complex_templates_df_filtered.filter(regex="Unnamed"), axis=1,
                                               inplace=True)
            complex_templates_df_filtered.reset_index(inplace=True, drop=True)

            if len(complex_templates_df_filtered) > self.predictor_config.max_template_count:
                break

        msa_out_path = outpath
        makedir_if_not_exists(msa_out_path)

        out_multimer_msas = []
        out_monomer_msas = []
        for chain_idx, chain_id in enumerate(chain_id_map):
            start_multimer_msa = os.path.join(start_msa_path, chain_id + '.start.multimer.a3m')
            fasta_chunks = (f">{chain_template_msas[chain_id]['desc'][i]}\n{chain_template_msas[chain_id]['seq'][i]}"
                            for i in range(len(chain_template_msas[chain_id]['desc'])))
            with open(start_multimer_msa + '.temp', 'w') as fw:
                fw.write('\n'.join(fasta_chunks) + '\n')
            
            out_multimer_msa = os.path.join(msa_out_path, f"{chain_id}.iteration{iteration}.multimer.a3m")
            if os.path.exists(start_multimer_msa):
                combine_a3ms([start_multimer_msa, f"{start_multimer_msa}.temp"], out_multimer_msa)
            else:
                os.system(f"cp {start_multimer_msa}.temp {out_multimer_msa}")
                
            out_multimer_msas += [out_multimer_msa]

            start_monomer_msa = os.path.join(start_msa_path, chain_id + ".start.monomer.a3m")
            if os.path.exists(start_monomer_msa):
                out_monomer_msa = os.path.join(msa_out_path, f"{chain_id}.iteration{iteration}.monomer.a3m")
                os.system("cp " + start_monomer_msa + " " + out_monomer_msa)
                out_monomer_msas += [out_monomer_msa]

        interact_dict = {}
        msa_len = -1
        for i in range(0, len(out_multimer_msas)):
            msa_sequences, msa_descriptions = parse_fasta(out_multimer_msas[i])
            current_len = len(msa_descriptions)
            if msa_len == -1:
                msa_len = current_len
            elif current_len != msa_len:
                raise Exception(f"The length of each msas are not equal! {out_multimer_msas}")
            interact_dict[f'index_{i + 1}'] = [j for j in range(msa_len)]
        interact_df = pd.DataFrame(interact_dict)
        interact_csv = os.path.join(outpath, f'interaction.iteration{iteration}.csv')
        interact_df.to_csv(interact_csv)

        top_template_files = []
        for template_result, chain_id in zip(template_results, chain_id_map):
            top_template_file = os.path.join(outpath, f"{chain_id}.top{self.predictor_config.max_template_count}")
            check_and_rank_monomer_templates_local_and_global(template_result=template_result,
                                                              outfile=top_template_file,
                                                              query_sequence=chain_id_map[chain_id].sequence,
                                                              max_template_count=self.predictor_config.max_template_count)
            top_template_files += [top_template_file]
        return top_template_files, out_monomer_msas, out_multimer_msas, interact_csv

    def copy_atoms_and_unzip(self, templates, outdir):
        os.chdir(outdir)
        num_templates = min(len(templates), 50)
        for i in range(num_templates):
            template_pdb = templates.loc[i, 'target']
            trg_pdb_path = os.path.join(outdir, template_pdb)
            print(trg_pdb_path)
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


    def search_single(self, chain_id_map, fasta_path, pdb_path, pkl_path, msa_paths, outdir):

        fasta_path = os.path.abspath(fasta_path)

        targetname = pathlib.Path(fasta_path).stem

        print(f"Processing {targetname}")

        outdir = os.path.abspath(outdir)

        makedir_if_not_exists(outdir)

        cwd = os.getcwd()

        ref_start_pdb = pdb_path
        ref_start_pkl = pkl_path
        ref_start_msa_paths = msa_paths
        ref_start_templates = {}

        model_iteration_scores = []

        print(f"Start to refine {ref_start_pdb}")

        common_parameters =   f"--fasta_path={fasta_path} " \
                              f"--env_dir={self.params['alphafold_env_dir']} " \
                              f"--database_dir={self.params['alphafold_database_dir']} " \
                              f"--benchmark={self.params['alphafold_benchmark']} " \
                              f"--use_gpu_relax={self.params['use_gpu_relax']} " \
                              f"--max_template_date={self.params['max_template_date']} "

        multimer_num_ensemble = self.get_heteromer_config(self.predictor_config, 'num_ensemble')
        multimer_num_recycle = self.get_heteromer_config(self.predictor_config, 'num_recycle')
        num_multimer_predictions_per_model = self.get_heteromer_config(self.predictor_config, 'predictions_per_model')
        model_preset = self.get_heteromer_config(self.predictor_config, 'model_preset')
        relax_topn_predictions = self.get_heteromer_config(self.predictor_config, 'relax_topn_predictions')
        dropout = self.get_heteromer_config(self.predictor_config, 'dropout')
        dropout_structure_module = self.get_heteromer_config(self.predictor_config, 'dropout_structure_module')

        common_parameters +=  f"--multimer_num_ensemble={multimer_num_ensemble} " \
                                f"--multimer_num_recycle={multimer_num_recycle} " \
                                f"--num_multimer_predictions_per_model={num_multimer_predictions_per_model} " \
                                f"--model_preset={model_preset} " \
                                f"--relax_topn_predictions={relax_topn_predictions} " \
                                f"--models_to_relax=TOPN "

        for num_iteration in range(self.predictor_config.max_iteration):
            os.chdir(cwd)
            current_work_dir = os.path.join(outdir, f"iteration{num_iteration + 1}")
            makedir_if_not_exists(current_work_dir)

            start_pdb = os.path.join(current_work_dir, "start.pdb")
            start_pkl = os.path.join(current_work_dir, "start.pkl")
            start_msa_path = os.path.join(current_work_dir, "start_msas")
            if os.path.exists(start_msa_path):
                os.system(f"rm -rf {start_msa_path}")
            makedir_if_not_exists(start_msa_path)

            if os.path.exists(ref_start_pkl):
                with open(ref_start_pkl, 'rb') as f:
                    ref_avg_lddt = float(pickle.load(f)['ranking_confidence'])
                os.system(f"cp {ref_start_pkl} {start_pkl}")
            else:
                ref_avg_lddt = 0

            for chain_id in ref_start_msa_paths:
                if len(ref_start_msa_paths[chain_id]['paired_msa']) > 0:
                    os.system(f"cp {ref_start_msa_paths[chain_id]['paired_msa']} " +
                            os.path.join(start_msa_path, chain_id+".start.multimer.a3m"))
                    os.system(f"cp {ref_start_msa_paths[chain_id]['monomer_msa']} " + 
                            os.path.join(start_msa_path, chain_id+".start.monomer.a3m"))

            os.system(f"cp {ref_start_pdb} {start_pdb}")

            model_iteration_scores += [ref_avg_lddt]

            out_model_dir = os.path.join(current_work_dir, "alphafold")

            if not complete_result(out_model_dir, 5 * int(num_multimer_predictions_per_model)):

                chain_pdbs = split_pdb_unrelax2relax(start_pdb, current_work_dir)

                template_results = []

                out_template_dir = os.path.join(current_work_dir, "templates")
                makedir_if_not_exists(out_template_dir)

                for chain_id in chain_pdbs:
                    print(chain_id)
                    if chain_id not in chain_id_map:
                        raise Exception("Multimer fasta file and model doesn't match!")

                    monomer_work_dir = os.path.join(current_work_dir, chain_id)
                    makedir_if_not_exists(monomer_work_dir)

                    chain_pdb = os.path.join(monomer_work_dir, chain_id+".pdb")
                    os.system(f"mv {chain_pdbs[chain_id]} {chain_pdb}")
                    foldseek_res = self.search_templates_foldseek(inpdb=chain_pdb, outdir=os.path.join(monomer_work_dir, 'foldseek'))

                    if len(foldseek_res['all_alignment']) == 0:
                        print(f"Cannot find any templates for {chain_id} in iteration {num_iteration + 1}")
                        break

                    template_results += [foldseek_res]

                    #self.copy_atoms_and_unzip(templates=foldseek_res['all_alignment'],
                    #                          outdir=out_template_dir)

                if len(template_results) != len(chain_id_map):
                    break

                template_files, monomer_msa_files, multimer_msa_files, msa_pair_file = \
                    self.concatenate_msa_and_templates(
                        chain_id_map=chain_id_map,
                        template_results=template_results,
                        start_msa_path=start_msa_path,
                        outpath=current_work_dir,
                        iteration=num_iteration + 1)

                find_templates = True
                for chain_id, template_file in zip(chain_id_map, template_files):
                    templates = pd.read_csv(template_file, sep='\t')
                    if len(templates) == 0:
                        print(f"Cannot find any templates for {chain_id} in iteration {num_iteration + 1}")
                        find_templates = False
                        continue
                    self.copy_atoms_and_unzip(templates=templates, outdir=out_template_dir)

                if not find_templates:
                    break

                makedir_if_not_exists(out_model_dir)

                cmd = f"python {self.params['alphafold_multimer_program']}  " \
                      f"--multimer_a3ms={','.join(multimer_msa_files)} " \
                      f"--msa_pair_file={msa_pair_file} " \
                      f"--monomer_temp_csvs={','.join(template_files)} " \
                      f"--struct_atom_dir={out_template_dir} " \
                      f"--output_dir={out_model_dir} " + common_parameters

                if len(monomer_msa_files) > 0:
                    cmd += f"--monomer_a3ms={','.join(monomer_msa_files)} "
                else:
                    cmd += f"--monomer_a3ms={','.join(multimer_msa_files)} "

                try:
                    os.chdir(self.params['alphafold_program_dir'])
                    print(cmd)
                    #os.system(cmd)
                except Exception as e:
                    print(e)

            new_ranking_json_file = os.path.join(out_model_dir, "ranking_debug.json")
            new_ranking_json = json.loads(open(new_ranking_json_file).read())
            max_lddt_score = new_ranking_json["iptm+ptm"][list(new_ranking_json["order"])[0]]

            print(f'#########Iteration: {num_iteration + 1}#############')
            print(f"plddt before: {ref_avg_lddt}")
            print(f"plddt after: {max_lddt_score}")
            if max_lddt_score > ref_avg_lddt:
                print("Continue to refine")
                ref_start_pdb = os.path.join(out_model_dir, "ranked_0.pdb")
                model_name = list(new_ranking_json["order"])[0]
                ref_start_pkl = os.path.join(out_model_dir, f"result_{model_name}.pkl") 
                ref_start_msa_paths = {}
                for chain_id in chain_id_map:
                    ref_start_msa_paths[chain_id] = dict(paired_msa=os.path.join(out_model_dir, f"msas/{chain_id}.paired.a3m"),
                                                         monomer_msa=os.path.join(out_model_dir, "msas", chain_id, "monomer_final.a3m"))
                    ref_start_templates[chain_id] = os.path.join(current_work_dir, f"{chain_id}.top{self.predictor_config.max_template_count}")

                print('##################################################')
                if num_iteration + 1 >= self.predictor_config.max_iteration:
                    print("Reach maximum iteration")
                    model_iteration_scores += [max_lddt_score]
            else:
                # keep the models in iteration 1 even through the plddt score decreases
                if num_iteration == 0:
                    ref_start_pdb = os.path.join(out_model_dir, "ranked_0.pdb")
                    model_name = list(new_ranking_json["order"])[0]
                    ref_start_pkl = os.path.join(out_model_dir, f"result_{model_name}.pkl") 
                    ref_start_msa_paths = {}
                    for chain_id in chain_id_map:
                        ref_start_msa_paths[chain_id] = dict(paired_msa=os.path.join(out_model_dir, f"msas/{chain_id}.paired.a3m"),
                                                         monomer_msa=os.path.join(out_model_dir, "msas", chain_id, "monomer_final.a3m"))
                        ref_start_templates[chain_id] = os.path.join(current_work_dir, f"{chain_id}.top{self.predictor_config.max_template_count}")
                    model_iteration_scores += [max_lddt_score]
                break

        # model_iteration_scores += [max_lddt_score]
        while len(model_iteration_scores) <= self.predictor_config.max_iteration:
            model_iteration_scores += [0]

        print(model_iteration_scores)
        df = pd.DataFrame(model_iteration_scores)
        df.to_csv(os.path.join(outdir, 'summary.csv'))

        final_model_dir = os.path.join(outdir, 'final')
        makedir_if_not_exists(final_model_dir)
        os.system(f"cp {ref_start_pdb} " + os.path.join(final_model_dir, "final.pdb"))
        os.system(f"cp {ref_start_pkl} " + os.path.join(final_model_dir, "final.pkl"))

        for chain_id in chain_id_map:
            os.system(f"cp {ref_start_msa_paths[chain_id]['paired_msa']} " +
                      os.path.join(final_model_dir, chain_id + ".paired.a3m"))
            os.system(f"cp {ref_start_msa_paths[chain_id]['monomer_msa']} " +
                      os.path.join(final_model_dir, chain_id + ".monomer.a3m"))
            os.system(f"cp {ref_start_templates[chain_id]} " +
                      os.path.join(final_model_dir, chain_id + ".templates"))

        os.chdir(cwd)

        return final_model_dir
