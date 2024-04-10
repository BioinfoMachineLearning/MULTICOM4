import copy
import os
import sys
import time
from multicom4.common.util import makedir_if_not_exists, check_dirs
import pandas as pd
from multiprocessing import Pool
import pathlib
from multicom4.common.protein import complete_result
import yaml, json
import numpy as np
from multicom4.common import config
from multicom4.monomer_structure_refinement import iterative_refine_pipeline
from multicom4.monomer_templates_search.structure_tmsearch_pipeline import monomer_tmsearch_based_template_search_pipeline

class Monomer_structure_prediction_pipeline_v2(config.pipeline):

    def __init__(self, params, run_methods=None):
        
        super().__init__()

        self.params = params

        self.non_af2_methods = ['paddle-helix', 'esmfold', 'deepfold', 'megafold']

        self.post_af2_methods = ['foldseek_refine', 'foldseek_refine_esm', 'foldseek_refine_esm_h']

        if run_methods is None:
            self.run_methods = ['default', 'default_seq_temp','def_drop_s','def_drop_nos',
                                'def_notemp', 'def_notemp_drop_s', 'def_notemp_drop_nos',
                                'original', 'ori_seq_temp', 'colabfold', 'colab_seq_temp',
                                'img', 'img_seq_temp', 'dhr', 
                                'deepmsa_dMSA_hhb', 'deepmsa_dMSA_jac', 'deepmsa_dMSA_hms',
                                'deepmsa_dMSA', 'deepmsa_qMSA', 'deepmsa_aMSA', 'deepmsa_qMSA_hhb',
                                'deepmsa_qMSA_jac', 'deepmsa_qMSA_hh3', 'deepmsa_qMSA_hms',
                                'deepmsa_DeepJGI_hms', 'deepmsa_DeepJGI', 'deepmsa_q3JGI', 
                                'deepmsa_q4JGI', 'deepmsa_q3JGI_hms', 'deepmsa_q4JGI_hms',
                                'paddle-helix', 'esmfold', 'deepfold', 'megafold',
                                'default_tmsearch', 'def_esm_msa',
                                'foldseek_refine', 'foldseek_refine_esm', 'foldseek_refine_esm_h']
        
        else:
            self.run_methods = run_methods

    def process_single(self, fasta_path, alndir, img_alndir, outdir, template_dir=None, run_script=False):
        
        targetname = pathlib.Path(fasta_path).stem

        cmd = ""

        makedir_if_not_exists(outdir)

        predictor_commands = {}
        
        for run_method in self.run_methods:
            
            cmds = []

            method_out_dir = os.path.join(outdir, run_method) 

            if run_method not in self.non_af2_methods and run_method not in self.post_af2_methods:
                
                predictor_config = self.monomer_config.predictors[run_method]
    
                monomer_num_ensemble = self.get_monomer_config(predictor_config, 'num_ensemble')
                monomer_num_recycle = self.get_monomer_config(predictor_config, 'num_recycle')
                num_monomer_predictions_per_model = self.get_monomer_config(predictor_config, 'predictions_per_model')
                model_preset = self.get_monomer_config(predictor_config, 'model_preset')
                relax_topn_predictions = self.get_monomer_config(predictor_config, 'relax_topn_predictions')
                dropout = self.get_monomer_config(predictor_config, 'dropout')
                dropout_structure_module = self.get_monomer_config(predictor_config, 'dropout_structure_module')

                common_parameters =  f"--fasta_path={fasta_path} " \
                                     f"--env_dir={self.params['alphafold_env_dir']} " \
                                     f"--database_dir={self.params['alphafold_database_dir']} " \
                                     f"--benchmark={self.params['alphafold_benchmark']} " \
                                     f"--use_gpu_relax={self.params['use_gpu_relax']} " \
                                     f"--max_template_date={self.params['max_template_date']} " \
                                     f"--monomer_num_ensemble={monomer_num_ensemble} " \
                                     f"--monomer_num_recycle={monomer_num_recycle} " \
                                     f"--num_monomer_predictions_per_model {num_monomer_predictions_per_model} " \
                                     f"--model_preset={model_preset} " \
                                     f"--relax_topn_predictions={relax_topn_predictions} " \
                                     f"--models_to_relax=TOPN "

                msa_source = self.get_monomer_config(predictor_config, 'msa_source')
                template_source = self.get_monomer_config(predictor_config, 'template_source')
                
                if msa_source == "default" and template_source == "pdb70":
                    errormsg = ""
                    bfd_uniref30_a3m = os.path.join(alndir, targetname + '_uniref30_bfd.a3m')
                    if not os.path.exists(bfd_uniref30_a3m):
                        errormsg = errormsg + f"Cannot find uniclust30 alignment for {targetname}: {bfd_uniref30_a3m}\n"
                    mgnify_sto = os.path.join(alndir, targetname + '_mgnify.sto')
                    if not os.path.exists(mgnify_sto):
                        errormsg = errormsg + f"Cannot find mgnify alignment for {targetname}: {mgnify_sto}\n"
                    uniref90_sto = os.path.join(alndir, targetname + '_uniref90.sto')
                    if not os.path.exists(uniref90_sto):
                        errormsg = errormsg + f"Cannot find uniref90 alignment for {targetname}: {uniref90_sto}\n"

                    if not complete_result(method_out_dir, 5 * num_monomer_predictions_per_model): 
                        if len(errormsg) > 0:
                            print(errormsg)
                        else:
                            cmd = f"cd {self.params['alphafold_program_dir']} && " \
                                  f"python {self.params['alphafold_default_program']} " \
                                  f"--bfd_uniref_a3ms={bfd_uniref30_a3m} " \
                                  f"--mgnify_stos={mgnify_sto} " \
                                  f"--uniref90_stos={uniref90_sto} " \
                                  f"--output_dir={method_out_dir} " + common_parameters
                            cmds += [cmd]
                
                else:

                    common_parameters += f"--dropout={dropout} " \
                                         f"--dropout_structure_module={dropout_structure_module} " \

                    errormsg = ""
                    if msa_source == "default":
                        bfd_uniref30_a3m = os.path.join(alndir, targetname + '_uniref30_bfd.a3m')
                        if not os.path.exists(bfd_uniref30_a3m):
                            errormsg = errormsg + f"Cannot find uniclust30 alignment for {targetname}: {bfd_uniref30_a3m}\n"
                        mgnify_sto = os.path.join(alndir, targetname + '_mgnify.sto')
                        if not os.path.exists(mgnify_sto):
                            errormsg = errormsg + f"Cannot find mgnify alignment for {targetname}: {mgnify_sto}\n"
                        uniref90_sto = os.path.join(alndir, targetname + '_uniref90.sto')
                        if not os.path.exists(uniref90_sto):
                            errormsg = errormsg + f"Cannot find uniref90 alignment for {targetname}: {uniref90_sto}\n"
                        common_parameters += f"--bfd_uniref_a3m={bfd_uniref30_a3m} " \
                                             f"--mgnify_sto={mgnify_sto} " \
                                             f"--uniref90_sto={uniref90_sto} " \

                    elif msa_source == "original":
                        uniref30_a3m = os.path.join(alndir, targetname + '_uniref30.a3m')
                        if not os.path.exists(uniref30_a3m):
                            errormsg = errormsg + f"Cannot find uniclust30 alignment for {targetname}: {uniref30_a3m}\n"
                        bfd_a3m = os.path.join(alndir, targetname + '_bfd.a3m')
                        if not os.path.exists(bfd_a3m):
                            errormsg = errormsg + f"Cannot find bfd alignment for {targetname}: {bfd_a3m}\n"
                        mgnify_sto = os.path.join(alndir, targetname + '_mgnify.sto')
                        if not os.path.exists(mgnify_sto):
                            errormsg = errormsg + f"Cannot find mgnify alignment for {targetname}: {mgnify_sto}\n"
                        uniref90_sto = os.path.join(alndir, targetname + '_uniref90.sto')
                        if not os.path.exists(uniref90_sto):
                            errormsg = errormsg + f"Cannot find uniref90 alignment for {targetname}: {uniref90_sto}\n"
                        common_parameters += f"--bfd_uniref_a3m={uniref30_a3m} " \
                                             f"--bfd_a3m={bfd_a3m} " \
                                             f"--mgnify_sto={mgnify_sto} " \
                                             f"--uniref90_sto={uniref90_sto} " \

                    elif msa_source == "colabfold":
                        colabfold_a3m = os.path.join(alndir, targetname + '_colabfold.a3m')
                        if not os.path.exists(colabfold_a3m):
                            errormsg = errormsg + f"Cannot find rosettafold alignment for {targetname}: {colabfold_a3m}\n"
                        common_parameters += f"--custom_msa={colabfold_a3m} "

                    elif msa_source == "img":
                        img_a3m = os.path.join(img_alndir, targetname + '.a3m')
                        if not os.path.exists(img_a3m):
                            errormsg = errormsg + f"Cannot find img alignment for {targetname}: {img_a3m}\n"
                        common_parameters += f"--custom_msa={img_a3m} "

                    elif msa_source == "dhr":
                        dhr_af_a3m = os.path.join(alndir, targetname + '_dhr_af.a3m')
                        if not os.path.exists(dhr_af_a3m):
                            errormsg = errormsg + f"Cannot find rosettafold alignment for {targetname}: {dhr_af_a3m}\n"
                        common_parameters += f"--custom_msa={dhr_af_a3m} "

                    elif msa_source == "esm_msa":

                        input_msa_path = os.path.join(outdir, predictor_config.input_msa_source, 'msas', 'monomer_final.a3m')

                        esm_msa_path = os.path.join(method_out_dir, 'esm.a3m')
                        
                        os.makedirs(method_out_dir, exist_ok=True)

                        if not os.path.exists(esm_msa_path):
                            cmd = f"sh {self.params['esm_msa_program']} {input_msa_path} {esm_msa_path}"
                            os.system(cmd)
                            if not os.path.exists(esm_msa_path):
                                print(f"Failed to generate the esm msa: {esm_msa_path}")
                        
                        common_parameters += f"--custom_msa={esm_msa_path} "

                    elif run_method.find('deepmsa_') >= 0:
                        deepmsa_a3m = os.path.join(alndir, 'DeepMSA2_a3m', 'finalMSAs', msa_source + '.a3m')
                        if not os.path.exists(deepmsa_a3m):
                            errormsg = errormsg + f"Cannot find deepmsa alignment for {targetname}: {deepmsa_a3m}\n"
                        common_parameters += f"--custom_msa={deepmsa_a3m} "
                    
                    if template_source == "pdb70":
                        uniref90_sto = os.path.join(alndir, targetname + '_uniref90.sto')
                        if not os.path.exists(uniref90_sto):
                            errormsg = errormsg + f"Cannot find uniref90 alignment for {targetname}: {uniref90_sto}\n"
                        if common_parameters.find('--uniref90_sto') < 0:
                            common_parameters += f"--uniref90_sto={uniref90_sto} "

                    elif template_source == "pdb_sort90":
                        temp_struct_csv = os.path.join(template_dir, "sequence_templates.csv")
                        if not os.path.exists(temp_struct_csv):
                            errormsg = errormsg + f"Cannot find template csv for {targetname}: {temp_struct_csv}\n"

                        struct_atom_dir = os.path.join(template_dir, "templates")
                        if not os.path.exists(struct_atom_dir):
                            errormsg = errormsg + f"Cannot find template directory for {targetname}: {struct_atom_dir}\n"
                        common_parameters += f"--temp_struct_csv={temp_struct_csv} " \
                                             f"--struct_atom_dir={struct_atom_dir} " \

                    elif template_source == "notemplate":
                        common_parameters += f"--notemplate=true "

                    elif template_source == "tmsearch":
                        default_ranked_0_pdb = os.path.join(outdir, 'default', 'ranked_0.pdb')
                        if not os.path.exists(default_ranked_0_pdb):
                            errormsg = errormsg + f"Cannot find ranked_0.pdb for default"
                        else:
                            workdir = os.path.join(method_out_dir, 'tmsearch')
                            makedir_if_not_exists(workdir)
                        
                            pipeline = monomer_tmsearch_based_template_search_pipeline(self.params)
                            template_file, template_dir = pipeline.search(fasta_path=fasta_path, inpdb=default_ranked_0_pdb, outdir=workdir)

                            common_parameters += f"--temp_struct_csv={template_file} " \
                                                 f"--struct_atom_dir={template_dir} " 

                    if not complete_result(method_out_dir, 5 * num_monomer_predictions_per_model): 
                        if len(errormsg) > 0:
                            print(errormsg)
                        else:
                            cmd = f"cd {self.params['alphafold_program_dir']} && " \
                                  f"python {self.params['alphafold_program']} " \
                                  f"--output_dir={method_out_dir} " + common_parameters
                            cmds += [cmd]

            elif run_method ==  "paddle-helix":
                os.makedirs(method_out_dir, exist_ok=True)
                if not os.path.exists(os.path.join(method_out_dir, 'unrelaxed.pdb')):
                    cmd = f"sh {self.params['paddle_helix_program']} {fasta_path} {method_out_dir}"
                    cmds += [cmd]

            elif run_method ==  "esmfold":
                os.makedirs(method_out_dir, exist_ok=True)
                for num_recycle in [4, 10, 50]:
                    outpdb = os.path.join(method_out_dir, f"{num_recycle}.pdb")
                    if not os.path.exists(outpdb):
                        cmd = f"sh {self.params['esmfold_program']} {fasta_path} {outpdb} {num_recycle}"
                        cmds += [cmd]

            elif run_method == "deepfold":
                os.makedirs(method_out_dir, exist_ok=True)
                pickle_path = os.path.join(outdir, 'default', 'features.pkl')
                # if not complete_result(method_out_dir, 5 * int(self.params['num_monomer_predictions_per_model'])):
                if not complete_result(method_out_dir, 5):
                    cmd = f"cd {self.params['deepfold_program_dir']} && " \
                          f"python {self.params['deepfold_program']} " \
                          f"--pickle_path {pickle_path} " \
                          f"--outdir {method_out_dir} "
                    cmds += [cmd]

            elif run_method == "megafold":
                os.makedirs(method_out_dir, exist_ok=True)
                a3m_path = os.path.join(alndir, 'colabfold_a3m', 'colabfold', 'search_result')
                out_a3m_path = os.path.join(method_out_dir, 'msas')
                if os.path.exists(out_a3m_path):
                    os.system(f"rm -rf {out_a3m_path}")
                os.makedirs(out_a3m_path)

                out_data_yaml = os.path.join(method_out_dir, 'megafold_data.yaml')
                yaml_template = self.params['megafold_data_yaml_template']
                with open(yaml_template, 'r') as f:
                    data = yaml.safe_load(f)

                data['database_search'] = {
                        'hhsearch_binary_path': os.path.join(self.params['alphafold_env_dir'], 'hhsearch'),
                        'kalign_binary_path': os.path.join(self.params['alphafold_env_dir'], 'kalign'),
                        'pdb70_database_path': self.params['pdb70_hhsuite_database'],
                        'mmcif_dir': os.path.join(self.params['alphafold_database_dir'], 'pdb_mmcif/mmcif_files'),
                        'obsolete_pdbs_path': os.path.join(self.params['alphafold_database_dir'], 'pdb_mmcif/obsolete.dat'),
                        'max_template_date': self.params['max_template_date'],
                        'mmseqs_binary': self.params['mmseq_program'],
                        'uniref30_path': os.path.join(self.params['colabfold_databases'], 'uniref30_2202_db'),
                        'database_envdb_dir': os.path.join(self.params['colabfold_databases'], 'colabfold_envdb_202108_db'),
                        'a3m_result_path': out_a3m_path,
                }

                with open(out_data_yaml, 'w') as fw:
                    yaml.safe_dump(data, fw)

                if not os.path.exists(os.path.join(method_out_dir, 'unrelaxed_seq.pdb')):
                    fastadir = os.path.join(method_out_dir, 'fasta')
                    os.makedirs(fastadir, exist_ok=True)
                    os.system(f"cp {fasta_path} {fastadir}/seq.fasta")
                    os.system(f"cp -r {a3m_path} {out_a3m_path}/seq")
                    os.system(f"rm {out_a3m_path}/seq/final.a3m {out_a3m_path}/seq/0.a3m")
                    cmd = f"sh {self.params['megafold_program']} {fastadir} {method_out_dir} {out_data_yaml}"
                    cmds += [cmd]     
                    

            elif run_method == "foldseek_refine" or run_method == "foldseek_refine_esm" or run_method == "foldseek_refine_esm_h":
                
                # refine default top-ranked models
                refinement_inputs = []
                default_workdir = os.path.join(outdir, 'default')
                ranking_json_file = os.path.join(default_workdir, "ranking_debug.json")
                if not os.path.exists(ranking_json_file):
                    continue
                ranking_json = json.loads(open(ranking_json_file).read())
                
                for i in range(self.monomer_config.predictors.foldseek_refine.number_of_input_models):
                    pdb_path = os.path.join(default_workdir, f"ranked_{i}.pdb")
                    model_name = list(ranking_json["order"])[i]
                    pkl_path = os.path.join(default_workdir, f"result_{model_name}.pkl")
                    msa_path = os.path.join(default_workdir, 'msas', "monomer_final.a3m")
                    refine_input = iterative_refine_pipeline.refinement_input(fasta_path=fasta_path,
                                                                              pdb_path=pdb_path,
                                                                              pkl_path=pkl_path,
                                                                              msa_path=msa_path)
                    refinement_inputs += [refine_input]

                refine_dir = os.path.join(method_out_dir, 'workdir')
                makedir_if_not_exists(refine_dir)
                pipeline = iterative_refine_pipeline.Monomer_iterative_refinement_pipeline_server(params=self.params, config_name=run_method)
                pipeline.search(refinement_inputs=refinement_inputs, outdir=refine_dir,
                                uniref90_sto=os.path.join(alndir, targetname + '_uniref90.sto'))

                final_dir = os.path.join(method_out_dir, 'finaldir')
                makedir_if_not_exists(final_dir)

                pipeline = iterative_refine_pipeline.Monomer_refinement_model_selection(params=self.params, config_name=run_method)
                pipeline.select_v1(indir=refine_dir, outdir=final_dir, prefix=run_method)
                pipeline.make_predictor_results(final_dir, method_out_dir)
            
            if len(cmds) > 0:
                predictor_commands[run_method] = cmds
       
        bash_files = []
        if os.path.exists(self.params['slurm_script_template']):
            bash_script_dir = os.path.join(outdir, 'slurm_scripts')
            os.makedirs(bash_script_dir, exist_ok=True)
            for predictor in predictor_commands:
                slurm_template_file = self.params['slurm_script_template']
                if predictor in self.non_af2_methods:
                    slurm_template_file = self.params['slurm_nonaf_script_template']

                bash_file = os.path.join(bash_script_dir, predictor + '.sh')
                print(f"Generating bash file for {predictor}: {bash_file}")
                name = open(fasta_path).readlines()[0].rstrip('\n')[1:]
                jobname = f"{name}_{predictor}"
                with open(bash_file, 'w') as fw:
                    for line in open(slurm_template_file):
                        line = line.replace("JOBNAME", jobname)
                        fw.write(line)
                    fw.write('\n'.join(predictor_commands[predictor]))
                bash_files += [bash_file]
            if run_script:
                for bash_file in bash_files:
                    os.system(f"sbatch {bash_file}")
        else:
            bash_script_dir = os.path.join(outdir, 'bash_scripts')
            os.makedirs(bash_script_dir, exist_ok=True)
            for predictor in predictor_commands:
                bash_file = os.path.join(bash_script_dir, predictor + '.sh')
                print(f"Generating bash file for {predictor}: {bash_file}")
                with open(bash_file, 'w') as fw:
                    fw.write('\n'.join(predictor_commands[predictor]))
                bash_files += [bash_file]
            
            if run_script:
                for bash_file in bash_files:
                    os.system(f"sh {bash_file}")

        # add ranking for deepmsa2 alignments
        deepmsa_alns = []
        deepmsa_plddts = []
        for method in self.run_methods:
            if method.find('deepmsa') >= 0:
                ranking_json_file = os.path.join(outdir, method, "ranking_debug.json")
                if not os.path.exists(ranking_json_file):
                    continue

                ranking_json = json.loads(open(ranking_json_file).read())
                plddts = []
                for key in ranking_json['plddts'].keys():
                    plddts += [float(ranking_json['plddts'][key])]
                deepmsa_alns += [method.replace('deepmsa_', '')]
                deepmsa_plddts += [np.max(np.array(plddts))]

        with open(os.path.join(outdir, 'deepmsa.rank'), 'w') as fw:
            indices = np.argsort(-np.array(deepmsa_plddts))
            for index in indices:
                fw.write(f"{deepmsa_alns[index]}\t{deepmsa_plddts[index]}\n")


    def process(self, monomers, alndir, outdir, templatedir=None):
        for fasta_path in monomers:
            fasta_name = pathlib.Path(fasta_path).stem
            monomer_aln_dir = os.path.join(alndir, fasta_name)
            monomer_outdir = os.path.join(outdir, fasta_name)
            monomer_template_dir = ""
            if templatedir is not None:
                monomer_template_dir = os.path.join(templatedir, fasta_name)
            self.process_single(fasta_path=fasta_path, alndir=monomer_aln_dir, outdir=monomer_outdir,
                                template_dir=monomer_template_dir)

        print("The tertiary structure generation for monomers has finished!")

    def process_single_domain(self, fasta_path, alndir, outdir, run_script=False):
        
        targetname = pathlib.Path(fasta_path).stem

        cmd = ""

        makedir_if_not_exists(outdir)

        predictor_commands = {}
        
        for run_method in self.run_methods:
            
            cmds = []

            method_out_dir = os.path.join(outdir, run_method) 

            if run_method not in self.non_af2_methods and run_method not in self.post_af2_methods:
                
                predictor_config = self.monomer_config.predictors[run_method]
    
                monomer_num_ensemble = self.get_monomer_config(predictor_config, 'num_ensemble')
                monomer_num_recycle = self.get_monomer_config(predictor_config, 'num_recycle')
                num_monomer_predictions_per_model = self.get_monomer_config(predictor_config, 'predictions_per_model')
                model_preset = self.get_monomer_config(predictor_config, 'model_preset')
                relax_topn_predictions = self.get_monomer_config(predictor_config, 'relax_topn_predictions')
                dropout = self.get_monomer_config(predictor_config, 'dropout')
                dropout_structure_module = self.get_monomer_config(predictor_config, 'dropout_structure_module')

                common_parameters =  f"--fasta_path={fasta_path} " \
                                     f"--env_dir={self.params['alphafold_env_dir']} " \
                                     f"--database_dir={self.params['alphafold_database_dir']} " \
                                     f"--benchmark={self.params['alphafold_benchmark']} " \
                                     f"--use_gpu_relax={self.params['use_gpu_relax']} " \
                                     f"--max_template_date={self.params['max_template_date']} " \
                                     f"--monomer_num_ensemble={monomer_num_ensemble} " \
                                     f"--monomer_num_recycle={monomer_num_recycle} " \
                                     f"--num_monomer_predictions_per_model {num_monomer_predictions_per_model} " \
                                     f"--model_preset={model_preset} " \
                                     f"--relax_topn_predictions={relax_topn_predictions} " \
                                     f"--models_to_relax=TOPN "

                msa_source = self.get_monomer_config(predictor_config, 'msa_source')
                template_source = self.get_monomer_config(predictor_config, 'template_source')
                
                if msa_source == "default" and template_source == "pdb70":
                    errormsg = ""
                    bfd_uniref30_a3m = os.path.join(alndir, targetname + '_uniref30_bfd.a3m')
                    if not os.path.exists(bfd_uniref30_a3m):
                        errormsg = errormsg + f"Cannot find uniclust30 alignment for {targetname}: {bfd_uniref30_a3m}\n"
                    mgnify_sto = os.path.join(alndir, targetname + '_mgnify.sto')
                    if not os.path.exists(mgnify_sto):
                        errormsg = errormsg + f"Cannot find mgnify alignment for {targetname}: {mgnify_sto}\n"
                    uniref90_sto = os.path.join(alndir, targetname + '_uniref90.sto')
                    if not os.path.exists(uniref90_sto):
                        errormsg = errormsg + f"Cannot find uniref90 alignment for {targetname}: {uniref90_sto}\n"

                    if not complete_result(method_out_dir, 5 * num_monomer_predictions_per_model): 
                        if len(errormsg) > 0:
                            print(errormsg)
                        else:
                            cmd = f"cd {self.params['alphafold_program_dir']} && " \
                                  f"python {self.params['alphafold_default_program']} " \
                                  f"--bfd_uniref_a3ms={bfd_uniref30_a3m} " \
                                  f"--mgnify_stos={mgnify_sto} " \
                                  f"--uniref90_stos={uniref90_sto} " \
                                  f"--output_dir={method_out_dir} " + common_parameters
                            cmds += [cmd]
                
                else:

                    common_parameters += f"--dropout={dropout} " \
                                         f"--dropout_structure_module={dropout_structure_module} " \

                    errormsg = ""
                    if msa_source == "default":
                        bfd_uniref30_a3m = os.path.join(alndir, targetname + '_uniref30_bfd.a3m')
                        if not os.path.exists(bfd_uniref30_a3m):
                            errormsg = errormsg + f"Cannot find uniclust30 alignment for {targetname}: {bfd_uniref30_a3m}\n"
                        mgnify_sto = os.path.join(alndir, targetname + '_mgnify.sto')
                        if not os.path.exists(mgnify_sto):
                            errormsg = errormsg + f"Cannot find mgnify alignment for {targetname}: {mgnify_sto}\n"
                        uniref90_sto = os.path.join(alndir, targetname + '_uniref90.sto')
                        if not os.path.exists(uniref90_sto):
                            errormsg = errormsg + f"Cannot find uniref90 alignment for {targetname}: {uniref90_sto}\n"
                        common_parameters += f"--bfd_uniref_a3m={bfd_uniref30_a3m} " \
                                             f"--mgnify_sto={mgnify_sto} " \
                                             f"--uniref90_sto={uniref90_sto} " \

                    elif msa_source == "original":
                        uniref30_a3m = os.path.join(alndir, targetname + '_uniref30.a3m')
                        if not os.path.exists(uniref30_a3m):
                            errormsg = errormsg + f"Cannot find uniclust30 alignment for {targetname}: {uniref30_a3m}\n"
                        bfd_a3m = os.path.join(alndir, targetname + '_bfd.a3m')
                        if not os.path.exists(bfd_a3m):
                            errormsg = errormsg + f"Cannot find bfd alignment for {targetname}: {bfd_a3m}\n"
                        mgnify_sto = os.path.join(alndir, targetname + '_mgnify.sto')
                        if not os.path.exists(mgnify_sto):
                            errormsg = errormsg + f"Cannot find mgnify alignment for {targetname}: {mgnify_sto}\n"
                        uniref90_sto = os.path.join(alndir, targetname + '_uniref90.sto')
                        if not os.path.exists(uniref90_sto):
                            errormsg = errormsg + f"Cannot find uniref90 alignment for {targetname}: {uniref90_sto}\n"
                        common_parameters += f"--bfd_uniref_a3m={uniref30_a3m} " \
                                             f"--bfd_a3m={bfd_a3m} " \
                                             f"--mgnify_sto={mgnify_sto} " \
                                             f"--uniref90_sto={uniref90_sto} " \

                    elif msa_source == "colabfold":
                        colabfold_a3m = os.path.join(alndir, targetname + '_colabfold.a3m')
                        if not os.path.exists(colabfold_a3m):
                            errormsg = errormsg + f"Cannot find colabfold alignment for {targetname}: {colabfold_a3m}\n"
                        common_parameters += f"--custom_msa={colabfold_a3m} "

                    elif msa_source == "img":
                        img_a3m = os.path.join(img_alndir, targetname + '.a3m')
                        if not os.path.exists(img_a3m):
                            errormsg = errormsg + f"Cannot find img alignment for {targetname}: {img_a3m}\n"
                        common_parameters += f"--custom_msa={img_a3m} "

                    elif msa_source == "dhr":
                        dhr_af_a3m = os.path.join(alndir, targetname + '_dhr_af.a3m')
                        if not os.path.exists(dhr_af_a3m):
                            errormsg = errormsg + f"Cannot find dhr alignment for {targetname}: {dhr_af_a3m}\n"
                        common_parameters += f"--custom_msa={dhr_af_a3m} "

                    elif msa_source == "esm_msa":

                        input_msa_path = os.path.join(outdir, predictor_config.input_msa_source, 'msas', 'monomer_final.a3m')

                        esm_msa_path = os.path.join(method_out_dir, 'esm.a3m')
                        
                        os.makedirs(method_out_dir, exist_ok=True)

                        if not os.path.exists(esm_msa_path):
                            cmd = f"sh {self.params['esm_msa_program']} {input_msa_path} {esm_msa_path}"
                            os.system(cmd)
                            if not os.path.exists(esm_msa_path):
                                errormsg = errormsg + f"Failed to generate the esm msa: {esm_msa_path}"
                        
                        common_parameters += f"--custom_msa={esm_msa_path} "

                    elif msa_source in ['dom_hhsearch', 'dom_parser', 'dom_unidoc', 'dom_manual']:
                        
                        dom_combine_a3m = os.path.join(method_outdir, 'domain_alns', 'final.a3m')

                        if not os.path.exists(dom_combine_a3m):

                            errormsg = errormsg + f"Cannot find domain alignment for {targetname}: {dom_combine_a3m}\n"

                        common_parameters += f"--custom_msa={dom_combine_a3m} "

                    elif run_method.find('deepmsa_') >= 0:
                        deepmsa_a3m = os.path.join(alndir, 'DeepMSA2_a3m', 'finalMSAs', msa_source + '.a3m')
                        if not os.path.exists(deepmsa_a3m):
                            errormsg = errormsg + f"Cannot find deepmsa alignment for {targetname}: {deepmsa_a3m}\n"
                        common_parameters += f"--custom_msa={deepmsa_a3m} "

                    if template_source == "pdb70":
                        uniref90_sto = os.path.join(alndir, targetname + '_uniref90.sto')
                        if not os.path.exists(uniref90_sto):
                            errormsg = errormsg + f"Cannot find uniref90 alignment for {targetname}: {uniref90_sto}\n"
                        if common_parameters.find('--uniref90_sto') < 0:
                            common_parameters += f"--uniref90_sto={uniref90_sto} "

                    elif template_source == "pdb_sort90":
                        temp_struct_csv = os.path.join(template_dir, "sequence_templates.csv")
                        if not os.path.exists(temp_struct_csv):
                            errormsg = errormsg + f"Cannot find template csv for {targetname}: {temp_struct_csv}\n"

                        struct_atom_dir = os.path.join(template_dir, "templates")
                        if not os.path.exists(struct_atom_dir):
                            errormsg = errormsg + f"Cannot find template directory for {targetname}: {struct_atom_dir}\n"
                        common_parameters += f"--temp_struct_csv={temp_struct_csv} " \
                                             f"--struct_atom_dir={struct_atom_dir} " \

                    elif template_source == "notemplate":
                        common_parameters += f"--notemplate=true "

                    elif template_source == "tmsearch":
                        default_ranked_0_pdb = os.path.join(outdir, 'default', 'ranked_0.pdb')
                        if not os.path.exists(default_ranked_0_pdb):
                            errormsg = errormsg + f"Cannot find ranked_0.pdb for default"
                        else:
                            workdir = os.path.join(method_out_dir, 'tmsearch')
                            makedir_if_not_exists(workdir)
                        
                            pipeline = monomer_tmsearch_based_template_search_pipeline(self.params)
                            template_file, template_dir = pipeline.search(fasta_path=fasta_path, inpdb=default_ranked_0_pdb, outdir=workdir)

                            common_parameters += f"--temp_struct_csv={template_file} " \
                                                 f"--struct_atom_dir={template_dir} " 

                    if not complete_result(method_out_dir, 5 * num_monomer_predictions_per_model): 
                        if len(errormsg) > 0:
                            print(errormsg)
                        else:
                            cmd = f"cd {self.params['alphafold_program_dir']} && " \
                                  f"python {self.params['alphafold_program']} " \
                                  f"--output_dir={method_out_dir} " + common_parameters
                            cmds += [cmd]

            elif run_method ==  "paddle-helix":
                os.makedirs(method_out_dir, exist_ok=True)
                if not os.path.exists(os.path.join(method_out_dir, 'unrelaxed.pdb')):
                    cmd = f"sh {self.params['paddle_helix_program']} {fasta_path} {method_out_dir}"
                    cmds += [cmd]

            elif run_method ==  "esmfold":
                os.makedirs(method_out_dir, exist_ok=True)
                for num_recycle in [4, 10, 50]:
                    outpdb = os.path.join(method_out_dir, f"{num_recycle}.pdb")
                    if not os.path.exists(outpdb):
                        cmd = f"sh {self.params['esmfold_program']} {fasta_path} {outpdb} {num_recycle}"
                        cmds += [cmd]

            elif run_method == "deepfold":
                os.makedirs(method_out_dir, exist_ok=True)
                pickle_path = os.path.join(outdir, 'default', 'features.pkl')
                # if not complete_result(method_out_dir, 5 * int(self.params['num_monomer_predictions_per_model'])):
                if not complete_result(method_out_dir, 5):
                    cmd = f"cd {self.params['deepfold_program_dir']} && " \
                          f"python {self.params['deepfold_program']} " \
                          f"--pickle_path {pickle_path} " \
                          f"--outdir {method_out_dir} "
                    cmds += [cmd]

            elif run_method == "megafold":
                os.makedirs(method_out_dir, exist_ok=True)
                a3m_path = os.path.join(alndir, 'colabfold_a3m', 'colabfold', 'search_result')
                out_a3m_path = os.path.join(method_out_dir, 'msas')
                if os.path.exists(out_a3m_path):
                    os.system(f"rm -rf {out_a3m_path}")
                os.makedirs(out_a3m_path)

                out_data_yaml = os.path.join(method_out_dir, 'megafold_data.yaml')
                yaml_template = self.params['megafold_data_yaml_template']
                with open(yaml_template, 'r') as f:
                    data = yaml.safe_load(f)

                data['database_search'] = {
                        'hhsearch_binary_path': os.path.join(self.params['alphafold_env_dir'], 'hhsearch'),
                        'kalign_binary_path': os.path.join(self.params['alphafold_env_dir'], 'kalign'),
                        'pdb70_database_path': self.params['pdb70_hhsuite_database'],
                        'mmcif_dir': os.path.join(self.params['alphafold_database_dir'], 'pdb_mmcif/mmcif_files'),
                        'obsolete_pdbs_path': os.path.join(self.params['alphafold_database_dir'], 'pdb_mmcif/obsolete.dat'),
                        'max_template_date': self.params['max_template_date'],
                        'mmseqs_binary': self.params['mmseq_program'],
                        'uniref30_path': os.path.join(self.params['colabfold_databases'], 'uniref30_2202_db'),
                        'database_envdb_dir': os.path.join(self.params['colabfold_databases'], 'colabfold_envdb_202108_db'),
                        'a3m_result_path': out_a3m_path,
                }

                with open(out_data_yaml, 'w') as fw:
                    yaml.safe_dump(data, fw)

                if not os.path.exists(os.path.join(method_out_dir, 'unrelaxed_seq.pdb')):
                    fastadir = os.path.join(method_out_dir, 'fasta')
                    os.makedirs(fastadir, exist_ok=True)
                    os.system(f"cp {fasta_path} {fastadir}/seq.fasta")
                    os.system(f"cp -r {a3m_path} {out_a3m_path}/seq")
                    os.system(f"rm {out_a3m_path}/seq/final.a3m {out_a3m_path}/seq/0.a3m")
                    cmd = f"sh {self.params['megafold_program']} {fastadir} {method_out_dir} {out_data_yaml}"
                    cmds += [cmd]     
            
            if len(cmds) > 0:
                predictor_commands[run_method] = cmds
       
        bash_files = []
        if os.path.exists(self.params['slurm_script_template']):
            bash_script_dir = os.path.join(outdir, 'slurm_scripts')
            os.makedirs(bash_script_dir, exist_ok=True)
            for predictor in predictor_commands:
                slurm_template_file = self.params['slurm_script_template']
                if predictor in self.non_af2_methods:
                    slurm_template_file = self.params['slurm_nonaf_script_template']

                bash_file = os.path.join(bash_script_dir, predictor + '.sh')
                print(f"Generating bash file for {predictor}: {bash_file}")
                name = open(fasta_path).readlines()[0].rstrip('\n')[1:]
                jobname = f"{name}_{predictor}"
                with open(bash_file, 'w') as fw:
                    for line in open(slurm_template_file):
                        line = line.replace("JOBNAME", jobname)
                        fw.write(line)
                    fw.write('\n'.join(predictor_commands[predictor]))
                bash_files += [bash_file]
            if run_script:
                for bash_file in bash_files:
                    os.system(f"sbatch {bash_file}")
        else:
            bash_script_dir = os.path.join(outdir, 'bash_scripts')
            os.makedirs(bash_script_dir, exist_ok=True)
            for predictor in predictor_commands:
                bash_file = os.path.join(bash_script_dir, predictor + '.sh')
                print(f"Generating bash file for {predictor}: {bash_file}")
                with open(bash_file, 'w') as fw:
                    fw.write('\n'.join(predictor_commands[predictor]))
                bash_files += [bash_file]
            
            if run_script:
                for bash_file in bash_files:
                    os.system(f"sh {bash_file}")

        # add ranking for deepmsa2 alignments
        deepmsa_alns = []
        deepmsa_plddts = []
        for method in self.run_methods:
            if method.find('deepmsa') >= 0:
                ranking_json_file = os.path.join(outdir, method, "ranking_debug.json")
                if not os.path.exists(ranking_json_file):
                    continue

                ranking_json = json.loads(open(ranking_json_file).read())
                plddts = []
                for key in ranking_json['plddts'].keys():
                    plddts += [float(ranking_json['plddts'][key])]
                deepmsa_alns += [method.replace('deepmsa_', '')]
                deepmsa_plddts += [np.max(np.array(plddts))]

        with open(os.path.join(outdir, 'deepmsa.rank'), 'w') as fw:
            indices = np.argsort(-np.array(deepmsa_plddts))
            for index in indices:
                fw.write(f"{deepmsa_alns[index]}\t{deepmsa_plddts[index]}\n")