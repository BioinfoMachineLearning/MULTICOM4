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


class Multimer_structure_prediction_pipeline_v2(config.pipeline):

    def __init__(self, params, run_methods=None):

        super().__init__()

        self.params = params

        self.non_af2_methods = ['esmfold']

        self.post_af2_methods = ['def_mul_refine']

        if run_methods is None:
            # non-alphafold2 predictors
            self.run_methods = ['default_multimer', 'def_mul_struct', 'def_mul_tmsearch', 'def_mul_pdb70',
                                'def_mul_pdb', 'def_mul_comp', 'def_mul_drop_s', #'def_mul_af',
                                'def_mul_drop_nos', 'def_mul_notemp', 'def_mul_not_drop_s', 
                                'def_mul_not_drop_nos', 'def_mul_nopair', 'def_mul_esm_msa', 'uniclust_ox_a3m',
                                'pdb_inter_ref_a3m', 'pdb_inter_ref_sto', 'pdb_inter_prot_sto', 
                                'unidist_ref_a3m', 'unidist_ref_sto', 'unidist_prot_sto', 'spec_inter_ref_a3m', 'spec_struct',
                                'spec_pdb70', 'spec_pdb', 'spec_comp', # 'spec_af',
                                'spec_inter_ref_sto', 'spec_inter_prot_sto', 
                                'str_inter_ref_a3m', 'str_inter_ref_sto', 'str_struct',
                                'str_pdb70', 'str_pdb', 'str_comp', 'str_inter_prot_sto', # 'str_af',
                                'AFProfile']

            self.run_methods += ['esmfold']
        else:
            self.run_methods = run_methods
            
    def process(self,
                fasta_path,
                chain_id_map,
                aln_dir,
                complex_aln_dir,
                template_dir,
                monomer_model_dir,
                output_dir,
                run_script=False):

        makedir_if_not_exists(output_dir)

        predictor_commands = {}

        for method in self.run_methods:
            
            cmds = []

            common_parameters = f"--fasta_path={fasta_path} " \
                                f"--env_dir={self.params['alphafold_env_dir']} " \
                                f"--database_dir={self.params['alphafold_database_dir']} " \
                                f"--benchmark={self.params['alphafold_benchmark']} " \
                                f"--use_gpu_relax={self.params['use_gpu_relax']} " \
                                f"--max_template_date={self.params['max_template_date']} "

            if method not in self.non_af2_methods and method not in self.post_af2_methods:
        
                outdir = os.path.join(output_dir, method)

                predictor_config = self.heteromer_config.predictors[method]

                multimer_num_ensemble = self.get_heteromer_config(predictor_config, 'num_ensemble')
                multimer_num_recycle = self.get_heteromer_config(predictor_config, 'num_recycle')
                num_multimer_predictions_per_model = self.get_heteromer_config(predictor_config, 'predictions_per_model')
                model_preset = self.get_heteromer_config(predictor_config, 'model_preset')
                relax_topn_predictions = self.get_heteromer_config(predictor_config, 'relax_topn_predictions')
                dropout = self.get_heteromer_config(predictor_config, 'dropout')
                dropout_structure_module = self.get_heteromer_config(predictor_config, 'dropout_structure_module')     

                common_parameters +=  f"--multimer_num_ensemble={multimer_num_ensemble} " \
                                        f"--multimer_num_recycle={multimer_num_recycle} " \
                                        f"--num_multimer_predictions_per_model={num_multimer_predictions_per_model} " \
                                        f"--model_preset={model_preset} " \
                                        f"--relax_topn_predictions={relax_topn_predictions} " \
                                        f"--models_to_relax=TOPN "

                msa_unpaired_source = self.get_heteromer_config(predictor_config, 'msa_unpaired_source')
                msa_paired_source = self.get_heteromer_config(predictor_config, 'msa_paired_source')
                template_source = self.get_heteromer_config(predictor_config, 'template_source')

                if method == "default_multimer":
                    # run alphafold default pipeline:
                    monomers = [chain_id for chain_id in chain_id_map]
                    if not complete_result(outdir, 5 * num_multimer_predictions_per_model):
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
                            cmd =  f"cd {self.params['alphafold_program_dir']} && " \
                                   f"python {self.params['alphafold_default_program']} " \
                                   f"--bfd_uniref_a3ms={','.join(bfd_uniref_a3ms)} " \
                                   f"--mgnify_stos={','.join(mgnify_stos)} " \
                                   f"--uniref90_stos={','.join(uniref90_stos)} " \
                                   f"--uniprot_stos={','.join(uniprot_stos)} " \
                                   f"--output_dir={outdir} " + common_parameters
                            cmds += [cmd]

                elif method == "AFProfile": 
                    
                    os.makedirs(outdir, exist_ok=True)

                    os.chdir(self.params['alphafold_program_dir'])

                    default_feature_pkl = os.path.join(output_dir, 'default_multimer', 'features.pkl')
                    os.system(f"cp {default_feature_pkl} {outdir}")
                    
                    predictor_config = self.heteromer_config.predictors[method]

                    confidence_threshold = self.get_heteromer_config(predictor_config, 'confidence_threshold')
                    max_iteration = self.get_heteromer_config(predictor_config, 'max_iteration')
                    learning_rate = self.get_heteromer_config(predictor_config, 'learning_rate')

                    cmd =  f"cd {self.params['alphafold_program_dir']} && " \
                           f"python {self.params['afprofile_program']} " \
                           f"--fasta_path={fasta_path} " \
                           f"--data_dir={self.params['alphafold_database_dir']} " \
                           f"--model_preset={model_preset} " \
                           f"--num_recycles={multimer_num_recycle} " \
                           f"--confidence_threshold={confidence_threshold} " \
                           f"--max_iter={max_iteration} " \
                           f"--learning_rate={learning_rate} " \
                           f"--output_dir={outdir} " \
                           f"--models_to_relax=TOPN " \
                           f"--relax_topn_predictions={relax_topn_predictions} "

                    if not complete_result(outdir, max_iteration):
                        cmds += [cmd]

                else:

                    common_parameters += f"--dropout={dropout} " \
                                            f"--dropout_structure_module={dropout_structure_module} " \

                    monomer_a3ms, multimer_a3ms = [], []

                    base_cmd =  f"cd {self.params['alphafold_program_dir']} && " \
                                f"python {self.params['alphafold_multimer_program']} "

                    # msa_unpaired_source
                    if msa_unpaired_source == "default":
                        for chain_id in chain_id_map:
                            monomer = chain_id
                            default_alphafold_monomer_a3m = os.path.join(output_dir, 'default_multimer', 'msas', chain_id, "monomer_final.a3m")
                            if not os.path.exists(default_alphafold_monomer_a3m):
                                raise Exception(f"Cannot find default alphafold alignments for {monomer}: {default_alphafold_monomer_a3m}")
                            monomer_a3ms += [default_alphafold_monomer_a3m]

                    base_cmd += f"--monomer_a3ms={','.join(monomer_a3ms)} "

                    # msa_paired_source
                    if msa_paired_source == "None":
                        base_cmd += "--msa_pairing_hetero=false "
                    else:
                        if msa_paired_source == "esm_msa":
                            msadir = os.path.join(outdir, 'esm_msas')
                            os.makedirs(msadir, exist_ok=True)

                            for chain_id in chain_id_map:
                                monomer = chain_id
                                default_alphafold_multimer_a3m = os.path.join(output_dir, 'default_multimer', 'msas', monomer + '.paired.a3m')
                                if not os.path.exists(default_alphafold_multimer_a3m):
                                    raise Exception(f"Cannot find default alphafold alignments for {monomer}: {default_alphafold_multimer_a3m}")
                                
                                esm_msa_path = os.path.join(msadir, f'{monomer}.esm.paired.a3m')
                                if not os.path.join(esm_msa_path):
                                    cmd = f"sh {self.params['esm_msa_program']} {default_alphafold_multimer_a3m} {esm_msa_path}"
                                    os.system(cmd)
                                    if not os.path.exists(esm_msa_path):
                                        print(f"Failed to generate the esm msa: {esm_msa_path}")

                                multimer_a3ms += [esm_msa_path]

                            msa_pair_file = os.path.join(msadir, "interact.csv")
                            interact_dict = {}
                            msa_len = -1
                            for i in range(len(multimer_a3ms)):
                                with open(multimer_a3ms[i]) as f:
                                    input_fasta_str = f.read()
                                msa_sequences, msa_descriptions = parse_fasta(input_fasta_str)
                                current_len = len(msa_descriptions)
                                if msa_len == -1:
                                    msa_len = current_len
                                elif current_len != msa_len:
                                    raise Exception(f"The length of each msas are not equal! {multimer_a3ms}")
                                interact_dict[f'index_{i + 1}'] = [j for j in range(msa_len)]
                            interact_df = pd.DataFrame(interact_dict)
                            interact_df.to_csv(msa_pair_file)

                        elif msa_paired_source == "default":
                            for chain_id in chain_id_map:
                                monomer = chain_id
                                default_alphafold_multimer_a3m = os.path.join(output_dir, 'default_multimer', 'msas', monomer + '.paired.a3m')
                                if not os.path.exists(default_alphafold_multimer_a3m):
                                    raise Exception(f"Cannot find default alphafold alignments for {monomer}: {default_alphafold_multimer_a3m}")
                                multimer_a3ms += [default_alphafold_multimer_a3m]

                            msa_pair_file = os.path.join(output_dir, "default_multimer/msas/interact.csv")
                            interact_dict = {}
                            msa_len = -1
                            for i in range(len(multimer_a3ms)):
                                with open(multimer_a3ms[i]) as f:
                                    input_fasta_str = f.read()
                                msa_sequences, msa_descriptions = parse_fasta(input_fasta_str)
                                current_len = len(msa_descriptions)
                                if msa_len == -1:
                                    msa_len = current_len
                                elif current_len != msa_len:
                                    raise Exception(f"The length of each msas are not equal! {multimer_a3ms}")
                                interact_dict[f'index_{i + 1}'] = [j for j in range(msa_len)]
                            interact_df = pd.DataFrame(interact_dict)
                            interact_df.to_csv(msa_pair_file)
                        else:
                            msa_pair_file = os.path.join(complex_aln_dir, msa_paired_source, msa_paired_source + "_interact.csv")
                            if len(pd.read_csv(msa_pair_file)) <= 1:
                                continue
                            multimer_a3ms = [os.path.join(complex_aln_dir, msa_paired_source, monomer + "_con.a3m") for monomer in monomers]

                        base_cmd += f"--multimer_a3ms={','.join(multimer_a3ms)} " \
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

                    elif template_source == "tmsearch_structure_based_template":
                        template_file = os.path.join(template_dir, "tmsearch", "structure_templates.csv")
                        struct_atom_dir = os.path.join(template_dir, "tmsearch", "templates")
                        base_cmd += f"--temp_struct_csv={template_file} " \
                                    f"--struct_atom_dir={struct_atom_dir} "

                    elif template_source == "sequence_based_template_pdb_sort90":
                        template_file = os.path.join(template_dir, "pdb_seq", "sequence_templates.csv")
                        struct_atom_dir = os.path.join(template_dir, "pdb_seq", "templates")
                        base_cmd += f"--temp_struct_csv={template_file} " \
                                    f"--struct_atom_dir={struct_atom_dir} "

                    elif template_source == "sequence_based_template_pdb_complex":
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

                    cmds += [base_cmd]
            
            elif method == "esmfold":
                # run esmfold
                method_out_dir = os.path.join(output_dir, "esmfold")
                os.makedirs(method_out_dir, exist_ok=True)
                for num_recycle in [4, 10, 50]:
                    outpdb = os.path.join(method_out_dir, f"{num_recycle}.pdb")
                    if not os.path.exists(outpdb):
                        cmd = f"sh {self.params['esmfold_program']} {fasta_path} {outpdb} {num_recycle}"
                        cmds += [cmd]

            elif method == "def_mul_refine":
                
                # refine default top-ranked models
                refinement_inputs = []
                default_workdir = os.path.join(outdir, 'default_multimer')
                ranking_json_file = os.path.join(default_workdir, "ranking_debug.json")
                if not os.path.exists(ranking_json_file):
                    continue
                ranking_json = json.loads(open(ranking_json_file).read())

                for i in range(self.get_heteromer_config.predictors.def_refine.number_of_input_models):
                    pdb_path = os.path.join(default_workdir, f"ranked_{i}.pdb")
                    model_name = list(ranking_json["order"])[i]
                    pkl_path = os.path.join(default_workdir, f"result_{model_name}.pkl")
                    msa_paths = {}
                    for chain_id in chain_id_map:
                        monomer_msa = os.path.join(default_workdir, 'msas', chain_id, "monomer_final.a3m")
                        paired_msa = os.path.join(default_workdir, 'msas', f"{chain_id}.paired.a3m")
                        msa_paths[chain_id] = dict(paired_msa=paired_msa,
                                                   monomer_msa=monomer_msa)

                    refine_input = iterative_refine_pipeline_multimer.refinement_input_multimer(chain_id_map=chain_id_map,
                                                                                                fasta_path=fasta_path,
                                                                                                pdb_path=pdb_path,
                                                                                                pkl_path=pkl_path,
                                                                                                msa_paths=msa_paths)
                    refinement_inputs += [refine_input]

                refine_dir = os.path.join(method_out_dir, 'workdir')
                makedir_if_not_exists(refine_dir)

                pipeline = iterative_refine_pipeline_multimer.Multimer_iterative_refinement_pipeline_server(params=params, config_name=run_method)
                pipeline.search(refinement_inputs=refinement_inputs, outdir=refine_dir, stoichiometry="heteromer")

                final_dir = os.path.join(method_out_dir, 'finaldir')
                makedir_if_not_exists(final_dir)

                pipeline = iterative_refine_pipeline_multimer.Multimer_refinement_model_selection()
                pipeline.select_v1(indir=refine_dir, outdir=final_dir)
                pipeline.make_predictor_results(final_dir, method_out_dir)


            predictor_commands[method] = cmds

        bash_script_dir = os.path.join(outdir, 'bash')
        os.makedirs(bash_script_dir, exist_ok=True)

        bash_files = []
        for predictor in predictor_commands:
            bash_file = os.path.join(bash_script_dir, predictor + '.sh')
            print(f"Generating bash file for {predictor}: {bash_file}")
            with open(bash_file, 'w') as fw:
                fw.write('\n'.join(predictor_commands[predictor]))
            bash_files += [bash_file]
        
        if run_script:
            for bash_file in bash_files:
                os.system(f"sh {bash_file}")

        print("The multimer structure generation for multimers has finished!")

    def process_deepmsa(self,
                        fasta_path,
                        chain_id_map,
                        aln_dir,
                        complex_aln_dir,
                        output_dir,
                        run_script=False):

        makedir_if_not_exists(output_dir)

        result_dirs = []

        common_parameters =   f"--fasta_path={fasta_path} " \
                              f"--env_dir={self.params['alphafold_env_dir']} " \
                              f"--database_dir={self.params['alphafold_database_dir']} " \
                              f"--benchmark={self.params['alphafold_benchmark']} " \
                              f"--use_gpu_relax={self.params['use_gpu_relax']} " \
                              f"--max_template_date={self.params['max_template_date']} "

        predictor_config = self.heteromer_config.predictors.deepmsa2
    
        multimer_num_ensemble = self.get_heteromer_config(predictor_config, 'num_ensemble')
        multimer_num_recycle = self.get_heteromer_config(predictor_config, 'num_recycle')
        num_multimer_predictions_per_model = self.get_heteromer_config(predictor_config, 'predictions_per_model')
        model_preset = self.get_heteromer_config(predictor_config, 'model_preset')
        relax_topn_predictions = self.get_heteromer_config(predictor_config, 'relax_topn_predictions')
        dropout = self.get_heteromer_config(predictor_config, 'dropout')
        dropout_structure_module = self.get_heteromer_config(predictor_config, 'dropout_structure_module')  
                
        common_parameters +=  f"--multimer_num_ensemble={multimer_num_ensemble} " \
                              f"--multimer_num_recycle={multimer_num_recycle} " \
                              f"--num_multimer_predictions_per_model {num_multimer_predictions_per_model} " \
                              f"--model_preset={model_preset} " \
                              f"--relax_topn_predictions={relax_topn_predictions} " \
                              f"--models_to_relax=TOPN "

        monomers = [chain_id for chain_id in chain_id_map]
            
        deepmsa_complex_aln_dir = os.path.join(complex_aln_dir, 'deepmsa2')

        ranking_file = os.path.join(deepmsa_complex_aln_dir, 'deepmsa_paired_ranking.csv')
        
        ranking_df = pd.read_csv(ranking_file)

        predictor_commands = {}

        for index, concatenate_method in enumerate(ranking_df.name):
            
            if concatenate_method.find('.csv') > 0:
                continue

            method_outdir = os.path.join(output_dir, f"deepmsa2_{index}")

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

            predictor_commands[f"deepmsa2_{index}"] = [base_cmd]

        bash_script_dir = os.path.join(outdir, 'bash')
        os.makedirs(bash_script_dir, exist_ok=True)

        bash_files = []
        for predictor in predictor_commands:
            bash_file = os.path.join(bash_script_dir, predictor + '.sh')
            print(f"Generating bash file for {predictor}: {bash_file}")
            with open(bash_file, 'w') as fw:
                fw.write('\n'.join(predictor_commands[predictor]))
            bash_files += [bash_file]
        
        if run_script:
            for bash_file in bash_files:
                os.system(f"sh {bash_file}")

        print("The multimer structure generation for multimers has finished!")