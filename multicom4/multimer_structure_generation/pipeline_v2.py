import copy
import os
import sys
import time
from multicom4.common.util import makedir_if_not_exists, check_dirs
from multicom4.common.protein import complete_result, parse_fasta
import pandas as pd
from multiprocessing import Pool
import pathlib
from multicom4.common import config

def get_complex_alignments_by_method(monomers, concatenate_method, aln_dir):
    a3ms_path = []
    for monomer in monomers:
        if concatenate_method == 'uniclust_oxmatch_a3m':
            monomer_a3m = os.path.join(aln_dir, monomer, f"{monomer}_uniclust30.a3m")
            if not os.path.exists(monomer_a3m):
                raise Exception(f"Cannot find alignment for {monomer}: {monomer_a3m}")
            a3ms_path += [monomer_a3m]
        elif concatenate_method.find('_uniref_a3m') > 0:
            monomer_a3m = os.path.join(aln_dir, monomer, f"{monomer}_uniref30.a3m")
            if not os.path.exists(monomer_a3m):
                raise Exception(f"Cannot find alignment for {monomer}: {monomer_a3m}")
            a3ms_path += [monomer_a3m]
        elif concatenate_method.find('_uniref_sto') > 0:
            monomer_a3m = os.path.join(aln_dir, monomer, f"{monomer}_uniref90.sto")
            if not os.path.exists(monomer_a3m):
                raise Exception(f"Cannot find alignment for {monomer}: {monomer_a3m}")
            a3ms_path += [monomer_a3m]
        elif concatenate_method.find('_uniprot_sto') > 0:
            monomer_a3m = os.path.join(aln_dir, monomer, f"{monomer}_uniprot.sto")
            if not os.path.exists(monomer_a3m):
                raise Exception(f"Cannot find alignment for {monomer}: {monomer_a3m}")
            a3ms_path += [monomer_a3m]
    return a3ms_path


class Multimer_structure_prediction_pipeline_v2:

    def __init__(self, params, run_methods=None):

        self.params = params
        self.configs = config.HETEROMULTIMER_CONFIG
        self.non_af2_methods = ['esmfold']

        if run_methods is None:
            # non-alphafold2 predictors
            self.run_methods = list(self.configs.predictors.keys())
            self.run_methods += ['esmfold']
        else:
            self.run_methods = run_methods

    def get_config(self, predictor_config, config_name):
        if config_name in predictor_config:
            return predictor_config[config_name]
        else:
            return self.configs.common_config[config_name]
            
    def process(self,
                fasta_path,
                chain_id_map,
                aln_dir,
                complex_aln_dir,
                template_dir,
                monomer_model_dir,
                output_dir):

        makedir_if_not_exists(output_dir)

        result_dirs, check_num_prediction_per_model = [], []

        common_parameters =   f"--fasta_path={fasta_path} " \
                              f"--env_dir={self.params['alphafold_env_dir']} " \
                              f"--database_dir={self.params['alphafold_database_dir']} " \
                              f"--benchmark={self.params['alphafold_benchmark']} " \
                              f"--use_gpu_relax={self.params['use_gpu_relax']} " \
                              f"--max_template_date={self.params['max_template_date']} "

        while True:
            # run esmfold
            method_out_dir = os.path.join(output_dir, "esmfold")
            os.makedirs(method_out_dir, exist_ok=True)
            for num_recycle in [4, 10, 50]:
                outpdb = os.path.join(method_out_dir, f"{num_recycle}.pdb")
                if not os.path.exists(outpdb):
                    cmd = f"sh {self.params['esmfold_program']} {fasta_path} {outpdb} {num_recycle}"
                    print(cmd)
                    os.system(cmd)

            for method in self.run_methods:
                
                outdir = os.path.join(output_dir, method)

                predictor_config = self.configs.predictors[run_method]
    
                multimer_num_ensemble = self.get_config(predictor_config, 'num_ensemble')
                multimer_num_recycle = self.get_config(predictor_config, 'num_recycle')
                num_multimer_predictions_per_model = self.get_config(predictor_config, 'predictions_per_model')
                model_preset = self.get_config(predictor_config, 'model_preset')
                relax_topn_predictions = self.get_config(predictor_config, 'relax_topn_predictions')
                dropout = self.get_config(predictor_config, 'dropout')
                dropout_structure_module = self.get_config(predictor_config, 'dropout_structure_module')     

                common_parameters +=  f"--multimer_num_ensemble={multimer_num_ensemble} " \
                                      f"--multimer_num_recycle={multimer_num_recycle} " \
                                      f"--num_multimer_predictions_per_model {num_multimer_predictions_per_model} " \
                                      f"--model_preset={model_preset} " \
                                      f"--relax_topn_predictions={relax_topn_predictions} " \
                                      f"--models_to_relax=TOPN "

                msa_unpaired_source = self.get_config(predictor_config, 'msa_unpaired_source')
                msa_paired_source = self.get_config(predictor_config, 'msa_paired_source')
                template_source = self.get_config(predictor_config, 'template_source')

                if run_method == "default_multimer":
                    # run alphafold default pipeline:
                    monomers = [chain_id for chain_id in chain_id_map]
                    if not complete_result(outdir, 5 * int(self.params['num_multimer_predictions_per_model'])):
                        os.chdir(self.params['alphafold_program_dir'])
                        bfd_uniref_a3ms = []
                        mgnify_stos = []
                        uniref90_stos = []
                        uniprot_stos = []
                        for chain_id in chain_id_map:
                            monomer = chain_id
                            monomer_bfd_uniref_a3m = os.path.join(aln_dir, monomer, f"{monomer}_uniref30_bfd.a3m")
                            if not os.path.exists(monomer_bfd_uniref_a3m):
                                raise Exception(f"Cannot find bfd and uniref30 a3m for {monomer}: {monomer_bfd_uniref_a3m}")
                            bfd_uniref_a3ms += [monomer_bfd_uniref_a3m]

                            monomer_mgnify_sto = os.path.join(aln_dir, monomer, f"{monomer}_mgnify.sto")
                            if not os.path.exists(monomer_mgnify_sto):
                                raise Exception(f"Cannot find mgnify sto for {monomer}: {monomer_mgnify_sto}")
                            mgnify_stos += [monomer_mgnify_sto]

                            monomer_uniref90_sto = os.path.join(aln_dir, monomer, f"{monomer}_uniref90.sto")
                            if not os.path.exists(monomer_uniref90_sto):
                                raise Exception(f"Cannot find uniref90 sto for {monomer}: {monomer_uniref90_sto}")
                            uniref90_stos += [monomer_uniref90_sto]

                            monomer_uniprot_sto = os.path.join(aln_dir, monomer, f"{monomer}_uniprot.sto")
                            if not os.path.exists(monomer_uniprot_sto):
                                raise Exception(f"Cannot find uniprot sto for {monomer}: {monomer_uniprot_sto}")
                            uniprot_stos += [monomer_uniprot_sto]

                        if not complete_result(outdir, 5 * num_multimer_predictions_per_model):
                            cmd = f"python {self.params['alphafold_default_program']} " \
                                  f"--bfd_uniref_a3ms={','.join(bfd_uniref_a3ms)} " \
                                  f"--mgnify_stos={','.join(mgnify_stos)} " \
                                  f"--uniref90_stos={','.join(uniref90_stos)} " \
                                  f"--uniprot_stos={','.join(uniprot_stos)} " \
                                  f"--output_dir={outdir} " + common_parameters

                            print(cmd)
                            os.system(cmd)

                    result_dirs += [outdir]
                    check_num_prediction_per_model += [num_multimer_predictions_per_model]

                else:

                    common_parameters += f"--dropout={dropout} " \
                                         f"--dropout_structure_module={dropout_structure_module} " \

                    monomer_a3ms, multimer_a3ms = [], []

                    # msa_unpaired_source
                    if msa_unpaired_source == "default":
                        for chain_id in chain_id_map:
                            monomer = chain_id
                            default_alphafold_monomer_a3m = os.path.join(output_dir, 'default_multimer', 'msas', chain_id, "monomer_final.a3m")
                            if not os.path.exists(default_alphafold_monomer_a3m):
                                raise Exception(f"Cannot find default alphafold alignments for {monomer}: {default_alphafold_monomer_a3m}")
                            monomer_a3ms += [default_alphafold_monomer_a3m]

                    # msa_paired_source
                    if msa_paired_source == "default":
                        for chain_id in chain_id_map:
                            monomer = chain_id
                            default_alphafold_multimer_a3m = os.path.join(output_dir, 'default_multimer', 'msas', monomer + '.paired.a3m')
                            if not os.path.exists(default_alphafold_multimer_a3m):
                                raise Exception(f"Cannot find default alphafold alignments for {monomer}: {default_alphafold_multimer_a3m}")
                            multimer_a3ms += [default_alphafold_multimer_a3m]

                        msa_pair_file = os.path.join(output_dir, "default_multimer/msas/interact.csv")
                        interact_dict = {}
                        msa_len = -1
                        for i in range(len(default_alphafold_multimer_a3ms)):
                            with open(default_alphafold_multimer_a3ms[i]) as f:
                                input_fasta_str = f.read()
                            msa_sequences, msa_descriptions = parse_fasta(input_fasta_str)
                            current_len = len(msa_descriptions)
                            if msa_len == -1:
                                msa_len = current_len
                            elif current_len != msa_len:
                                raise Exception(f"The length of each msas are not equal! {default_alphafold_multimer_a3ms}")
                            interact_dict[f'index_{i + 1}'] = [j for j in range(msa_len)]
                        interact_df = pd.DataFrame(interact_dict)
                        interact_df.to_csv(msa_pair_file)

                    else:
                        msa_pair_file = os.path.join(complex_aln_dir, concatenate_method, concatenate_method + "_interact.csv")
                        if len(pd.read_csv(msa_pair_file)) <= 1:
                            continue
                        multimer_a3ms = [os.path.join(complex_aln_dir, concatenate_method, monomer + "_con.a3m") for monomer in monomers]

                    base_cmd = f"python {self.params['alphafold_multimer_program']} " \
                               f"--monomer_a3ms={','.join(monomer_a3ms)} " \
                               f"--multimer_a3ms={','.join(multimer_a3ms)} " \
                               f"--msa_pair_file={msa_pair_file} " \

                    # template_source
                    if template_source == "pdb_seqres":
                        template_stos = []
                        for chain_id in chain_id_map:
                            monomer = chain_id
                            monomer_template_sto = os.path.join(aln_dir, monomer, f"{monomer}_uniref90.sto")
                            if not os.path.exists(monomer_template_sto):
                                raise Exception(f"Cannot find template stos for {monomer}: {monomer_template_sto}")
                            template_stos += [monomer_template_sto]
                        base_cmd += f"--template_stos {','.join(template_stos)} "

                    elif template_source == "foldseek_structure_based_template":
                        template_file = os.path.join(template_dir, "struct_temp", "structure_templates.csv")
                        struct_atom_dir = os.path.join(template_dir, "struct_temp", "templates")
                        base_cmd += f"--temp_struct_csv={template_file} " \
                                    f"--struct_atom_dir={struct_atom_dir} "

                    elif template_source == "sequence_based_template_pdb":
                        template_file = os.path.join(template_dir, "pdb_seq", "sequence_templates.csv")
                        struct_atom_dir = os.path.join(template_dir, "pdb_seq", "templates")
                        base_cmd += f"--temp_struct_csv={template_file} " \
                                    f"--struct_atom_dir={struct_atom_dir} "

                    elif template_source == "sequence_based_template_complex_pdb":
                        template_file = os.path.join(template_dir, "complex_pdb_seq", "sequence_templates.csv")
                        struct_atom_dir = os.path.join(template_dir, "complex_pdb_seq", "templates")
                        base_cmd += f"--temp_struct_csv={template_file} " \
                                    f"--struct_atom_dir={struct_atom_dir} "

                    elif template_source == "sequence_based_template_pdb70":
                        template_file = os.path.join(template_dir, "pdb70_seq", "sequence_templates.csv")
                        template_hits_files = []
                        for monomer in monomers:
                            template_hits_file = os.path.join(template_dir, "pdb70_seq", monomer, "output.hhr")
                            if not os.path.exists(template_hits_file):
                                raise Exception(f"Cannot find template hit file for {monomer}: {template_hits_file}")
                            template_hits_files += [template_hits_file]
                        base_cmd += f"--temp_seq_pair_file={template_file} " \
                                    f"--template_hits_files={','.join(template_hits_files)} "

                    elif template_source == "alphafold_model_templates":
                        monomer_paths = []
                        for monomer in monomers:
                            monomer_path = os.path.join(monomer_model_dir, monomer, "default")
                            if not os.path.exists(monomer_path):
                                raise Exception(f"Cannot find monomer directory for {monomer}: {monomer_path}")
                            monomer_paths += [monomer_path]
                        base_cmd += f"--monomer_model_paths={','.join(monomer_paths)} "

                    elif template_source == "notemplate":
                        base_cmd += "--notemplate=true "  

                    base_cmd += f"--output_dir={outdir} " + common_parameters

                    if complete_result(outdir, 5 * num_multimer_predictions_per_model):
                        continue

                    makedir_if_not_exists(outdir)

                    print(base_cmd)
                    os.system(base_cmd)

                    result_dirs += [outdir]
                    check_num_prediction_per_model += [num_multimer_predictions_per_model]

            rerun = False
            for result_dir, num_prediction_per_model in zip(result_dirs, check_num_prediction_per_model):
                if not complete_result(result_dir, 5 * num_prediction_per_model):
                    rerun = True

            if not rerun:
                break

            print("end")

        print("The multimer structure generation for multimers has finished!")

    def process_deepmsa(self,
                        fasta_path,
                        chain_id_map,
                        aln_dir,
                        complex_aln_dir,
                        output_dir):

        makedir_if_not_exists(output_dir)

        result_dirs = []

        common_parameters =   f"--fasta_path={fasta_path} " \
                              f"--env_dir={self.params['alphafold_env_dir']} " \
                              f"--database_dir={self.params['alphafold_database_dir']} " \
                              f"--benchmark={self.params['alphafold_benchmark']} " \
                              f"--use_gpu_relax={self.params['use_gpu_relax']} " \
                              f"--max_template_date={self.params['max_template_date']} "

        predictor_config = self.configs.predictors.deepmsa2
    
        multimer_num_ensemble = self.get_config(predictor_config, 'num_ensemble')
        multimer_num_recycle = self.get_config(predictor_config, 'num_recycle')
        num_multimer_predictions_per_model = self.get_config(predictor_config, 'predictions_per_model')
        model_preset = self.get_config(predictor_config, 'model_preset')
        relax_topn_predictions = self.get_config(predictor_config, 'relax_topn_predictions')
        dropout = self.get_config(predictor_config, 'dropout')
        dropout_structure_module = self.get_config(predictor_config, 'dropout_structure_module')  
                
        common_parameters +=  f"--multimer_num_ensemble={multimer_num_ensemble} " \
                              f"--multimer_num_recycle={multimer_num_recycle} " \
                              f"--num_multimer_predictions_per_model {num_multimer_predictions_per_model} " \
                              f"--model_preset={model_preset} " \
                              f"--relax_topn_predictions={relax_topn_predictions} " \
                              f"--models_to_relax=TOPN "

        monomers = [chain_id for chain_id in chain_id_map]

        while True:
            
            deepmsa_complex_aln_dir = os.path.join(complex_aln_dir, 'deepmsa_species')

            for concatenate_method in os.listdir(deepmsa_complex_aln_dir):
                
                if concatenate_method.find('.csv') > 0:
                    continue

                method_outdir = os.path.join(output_dir, concatenate_method)

                os.chdir(self.params['alphafold_program_dir'])

                msa_pair_file = os.path.join(deepmsa_complex_aln_dir, concatenate_method, concatenate_method + "_interact.csv")
                if len(pd.read_csv(msa_pair_file)) <= 1:
                    continue

                paired_a3m_paths = [os.path.join(deepmsa_complex_aln_dir, concatenate_method, monomer + "_con.a3m") for monomer in monomers]
                monomer_a3m_paths = []
                template_stos = []

                monomer_a3m_names = concatenate_method.split('_')
                for chain_idx, chain_id in enumerate(chain_id_map):

                    monomer = chain_id

                    monomer_a3m_path = os.path.join(aln_dir, monomer, 'DeepMSA2_a3m', 'finalMSAs', monomer_a3m_names[chain_idx] + '.a3m')
                    if not os.path.exists(monomer_a3m_path):
                        raise Exception(f"Cannot find deepmsa a3m for {monomer}: {monomer_a3m_path}")

                    monomer_a3m_paths += [monomer_a3m_path]

                    monomer_template_sto = os.path.join(aln_dir, monomer, f"{monomer}_uniref90.sto")
                    if not os.path.exists(monomer_template_sto):
                        raise Exception(f"Cannot find template stos for {monomer}: {monomer_template_sto}")
                    template_stos += [monomer_template_sto]

                base_cmd = f"python {self.params['alphafold_multimer_program']} " \
                           f"--monomer_a3ms={','.join(monomer_a3m_paths)} " \
                           f"--multimer_a3ms={','.join(paired_a3m_paths)} " \
                           f"--msa_pair_file={msa_pair_file} " \
                           f"--template_stos {','.join(template_stos)} " \
                           f"--output_dir={method_outdir} " + common_parameters

                if complete_result(method_outdir, 5 * num_multimer_predictions_per_model):
                    continue

                makedir_if_not_exists(method_outdir)

                print(base_cmd)
                os.system(base_cmd)

                result_dirs += [method_outdir]

            rerun = False
            for result_dir in result_dirs:
                if not complete_result(result_dir, 5 * num_multimer_predictions_per_model):
                    rerun = True

            if not rerun:
                break

            print("end")

        print("The multimer structure generation for multimers has finished!")