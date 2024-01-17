import copy
import os
import sys
import time
from multicom4.common.util import makedir_if_not_exists, check_dirs
import pandas as pd
from multiprocessing import Pool
import pathlib
from multicom4.common.protein import complete_result
import yaml

class Monomer_structure_prediction_pipeline_v2:

    def __init__(self, params, run_methods=None):

        self.params = params
        self.deepmsa_noimg_tags = ['dMSA.hhb', 'dMSA.jac', 'dMSA.hms', 'dMSA', 'qMSA', 'aMSA', 'qMSA.hhb', 'qMSA.jac', 'qMSA.hh3', 'qMSA.hms']
        self.deepmsa_img_tags = ['DeepJGI.hms', 'DeepJGI', 'q3JGI', 'q4JGI', 'q3JGI.hms', 'q4JGI.hms']

        if run_methods is None:
            #template/no_template, nodropout/dropout dropout->struct_dropout/no
            # default: default_template_nodropout_struct_dropout
            self.run_methods = ['default', 'default+seq_template', 
                                'default+template+dropout+struct_dropout',
                                'default+template+dropout+no_struct_dropout',
                                'default+notemplate+nodropout+struct_dropout',
                                'default+notemplate+dropout+struct_dropout',
                                'default+notemplate+dropout+no_struct_dropout',

                                'original', 'original+seq_template',
                                'colabfold', 'colabfold+seq_template',
                                'dhr',
                                'paddle-helix', 'esmfold',
                                'deepfold', 'megafold']
            for deepmsa_noimg_tag in self.deepmsa_noimg_tags:
                self.run_methods += ['deepmsa_noimg_' + deepmsa_noimg_tag]
            for deepmsa_img_tag in self.deepmsa_img_tags:
                self.run_methods += ['deepmsa_img_' + deepmsa_img_tag]

        else:
            self.run_methods = run_methods

        self.method2dir = {'default': 'default',
                           'default+seq_template': 'default_seq_temp',

                           'default+template+dropout+struct_dropout': 'def_drop_s',
                           'default+template+dropout+no_struct_dropout': 'def_drop_nos',
                           'default+notemplate+nodropout+struct_dropout': 'def_notemp',
                           'default+notemplate+dropout+struct_dropout': 'def_notemp_drop_s',
                           'default+notemplate+dropout+no_struct_dropout': 'def_notemp_drop_nos',

                           'original': 'original',
                           'original+seq_template': 'ori_seq_temp',
                           'colabfold': 'colabfold',
                           'colabfold+seq_template': 'colab_seq_temp',
                           'img': 'img',
                           'img_seq_template': 'img_seq_temp',
                           'dhr': 'dhr',
                           'paddle-helix': 'paddle-helix',
                           'esmfold': 'esmfold',
                           'deepfold': 'deepfold',
                           'megafold': 'megafold',
                           'deepmsa_noimg_dMSA.hhb': 'deepmsa_dMSA_hhb',
                           'deepmsa_noimg_dMSA.jac': 'deepmsa_dMSA_jac',
                           'deepmsa_noimg_dMSA.hms': 'deepmsa_dMSA_hms',
                           'deepmsa_noimg_dMSA': 'deepmsa_dMSA',
                           'deepmsa_noimg_qMSA': 'deepmsa_qMSA',
                           'deepmsa_noimg_aMSA': 'deepmsa_aMSA',
                           'deepmsa_noimg_qMSA.hhb': 'deepmsa_qMSA_hhb',
                           'deepmsa_noimg_qMSA.jac': 'deepmsa_qMSA_jac',
                           'deepmsa_noimg_qMSA.hh3': 'deepmsa_qMSA_hh3',
                           'deepmsa_noimg_qMSA.hms': 'deepmsa_qMSA_hms',
                           'deepmsa_img_DeepJGI.hms': 'deepmsa_DeepJGI_hms', 
                           'deepmsa_img_DeepJGI': 'deepmsa_DeepJGI',  
                           'deepmsa_img_q3JGI': 'deepmsa_q3JGI', 
                           'deepmsa_img_q4JGI': 'deepmsa_q4JGI', 
                           'deepmsa_img_q3JGI.hms': 'deepmsa_q3JGI_hms', 
                           'deepmsa_img_q4JGI.hms': 'deepmsa_q4JGI_hms',
        }


    def process_single(self, fasta_path, alndir, outdir, template_dir=None):
        
        targetname = pathlib.Path(fasta_path).stem

        cmd = ""

        makedir_if_not_exists(outdir)

        common_parameters =   f"--fasta_path={fasta_path} " \
                              f"--env_dir={self.params['alphafold_env_dir']} " \
                              f"--database_dir={self.params['alphafold_database_dir']} " \
                              f"--monomer_num_ensemble={self.params['monomer_num_ensemble']} " \
                              f"--monomer_num_recycle={self.params['monomer_num_recycle']} " \
                              f"--num_monomer_predictions_per_model {self.params['num_monomer_predictions_per_model']} " \
                              f"--model_preset={self.params['monomer_model_preset']} " \
                              f"--benchmark={self.params['alphafold_benchmark']} " \
                              f"--use_gpu_relax={self.params['use_gpu_relax']} " \
                              f"--models_to_relax={self.params['models_to_relax']} " \
                              f"--max_template_date={self.params['max_template_date']} "

        for run_method in self.run_methods:
            
            cmds = []

            method_out_dir = os.path.join(outdir, self.method2dir[run_method]) 

            if run_method == "default":
                
                os.chdir(self.params['alphafold_program_dir'])

                errormsg = ""

                if not os.path.exists(alndir):
                    errormsg = errormsg + f"Cannot find alignment directory for {targetname}: {alndir}\n"

                bfd_uniref30_a3m = os.path.join(alndir, targetname + '_uniref30_bfd.a3m')
                if not os.path.exists(bfd_uniref30_a3m):
                    errormsg = errormsg + f"Cannot find uniclust30 alignment for {targetname}: {bfd_uniref30_a3m}\n"

                mgnify_sto = os.path.join(alndir, targetname + '_mgnify.sto')
                if not os.path.exists(mgnify_sto):
                    errormsg = errormsg + f"Cannot find mgnify alignment for {targetname}: {mgnify_sto}\n"

                uniref90_sto = os.path.join(alndir, targetname + '_uniref90.sto')
                if not os.path.exists(uniref90_sto):
                    errormsg = errormsg + f"Cannot find uniref90 alignment for {targetname}: {uniref90_sto}\n"

                if len(errormsg) == 0:
                    if not complete_result(method_out_dir, 5 * int(self.params['num_monomer_predictions_per_model'])):
                        cmd = f"python {self.params['alphafold_default_program']} " \
                            f"--bfd_uniref_a3ms={bfd_uniref30_a3m} " \
                            f"--mgnify_stos={mgnify_sto} " \
                            f"--uniref90_stos={uniref90_sto} " \
                            f"--output_dir={method_out_dir} " + common_parameters
                        cmds += [cmd]
                else:
                    print(errormsg)

            elif run_method == "default+seq_template":
                
                os.chdir(self.params['alphafold_program_dir'])

                errormsg = ""

                if not os.path.exists(alndir):
                    errormsg = errormsg + f"Cannot find alignment directory for {targetname}: {alndir}\n"

                bfd_uniref30_a3m = os.path.join(alndir, targetname + '_uniref30_bfd.a3m')
                if not os.path.exists(bfd_uniref30_a3m):
                    errormsg = errormsg + f"Cannot find uniclust30 alignment for {targetname}: {bfd_uniref30_a3m}\n"

                mgnify_sto = os.path.join(alndir, targetname + '_mgnify.sto')
                if not os.path.exists(mgnify_sto):
                    errormsg = errormsg + f"Cannot find mgnify alignment for {targetname}: {mgnify_sto}\n"

                uniref90_sto = os.path.join(alndir, targetname + '_uniref90.sto')
                if not os.path.exists(uniref90_sto):
                    errormsg = errormsg + f"Cannot find uniref90 alignment for {targetname}: {uniref90_sto}\n"

                if len(errormsg) == 0:
                    if not complete_result(method_out_dir, 5 * int(self.params['num_monomer_predictions_per_model'])):
                        temp_struct_csv = os.path.join(template_dir, "sequence_templates.csv")
                        struct_atom_dir = os.path.join(template_dir, "templates")
                        cmd = f"python {self.params['alphafold_program']} " \
                            f"--bfd_uniref_a3m={bfd_uniref30_a3m} " \
                            f"--mgnify_sto={mgnify_sto} " \
                            f"--uniref90_sto={uniref90_sto} " \
                            f"--temp_struct_csv={temp_struct_csv} " \
                            f"--struct_atom_dir={struct_atom_dir} " \
                            f"--output_dir={method_out_dir} " + common_parameters
                        cmds += [cmd]
                else:
                    print(errormsg)

            elif run_method in ["default+template+dropout+struct_dropout",
                                "default+template+dropout+no_struct_dropout",
                                'default+notemplate+nodropout+struct_dropout',
                                'default+notemplate+dropout+struct_dropout',
                                'default+notemplate+dropout+no_struct_dropout']:

                os.chdir(self.params['alphafold_program_dir'])

                errormsg = ""

                if not os.path.exists(alndir):
                    errormsg = errormsg + f"Cannot find alignment directory for {targetname}: {alndir}\n"

                bfd_uniref30_a3m = os.path.join(alndir, targetname + '_uniref30_bfd.a3m')
                if not os.path.exists(bfd_uniref30_a3m):
                    errormsg = errormsg + f"Cannot find uniclust30 alignment for {targetname}: {bfd_uniref30_a3m}\n"

                mgnify_sto = os.path.join(alndir, targetname + '_mgnify.sto')
                if not os.path.exists(mgnify_sto):
                    errormsg = errormsg + f"Cannot find mgnify alignment for {targetname}: {mgnify_sto}\n"

                uniref90_sto = os.path.join(alndir, targetname + '_uniref90.sto')
                if not os.path.exists(uniref90_sto):
                    errormsg = errormsg + f"Cannot find uniref90 alignment for {targetname}: {uniref90_sto}\n"

                if len(errormsg) == 0:
                    configs = run_method.split('+')
                    notemplate=False
                    dropout=False
                    dropout_structure_module=True
                    if configs[1] == "notemplate":
                        notemplate = True
                    if configs[2] == "dropout":
                        dropout=True
                    if configs[3] == "no_struct_dropout":
                        dropout_structure_module = False

                    if not complete_result(method_out_dir, 5 * int(self.params['num_monomer_predictions_per_model'])):
                        cmd = f"python {self.params['alphafold_program']} " \
                            f"--bfd_uniref_a3m={bfd_uniref30_a3m} " \
                            f"--mgnify_sto={mgnify_sto} " \
                            f"--uniref90_sto={uniref90_sto} " \
                            f"--notemplate={notemplate} " \
                            f"--dropout={dropout} " \
                            f"--dropout_structure_module={dropout_structure_module} " \
                            f"--output_dir={method_out_dir} " + common_parameters
                        cmds += [cmd]
                else:
                    print(errormsg)
                    
            elif run_method == "original":

                os.chdir(self.params['alphafold_program_dir'])

                errormsg = ""

                if not os.path.exists(alndir):
                    errormsg = errormsg + f"Cannot find alignment directory for {targetname}: {alndir}\n"

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

                if len(errormsg) == 0:
                    if not complete_result(method_out_dir, 5 * int(self.params['num_monomer_predictions_per_model'])):
                        cmd = f"python {self.params['alphafold_program']} " \
                            f"--bfd_uniref_a3m={uniref30_a3m} " \
                            f"--bfd_a3m={bfd_a3m} " \
                            f"--mgnify_sto={mgnify_sto} " \
                            f"--uniref90_sto={uniref90_sto} " \
                            f"--output_dir={method_out_dir} " + common_parameters
                        cmds += [cmd]
                else:
                    print(errormsg)

            elif run_method ==  "original+seq_template":

                os.chdir(self.params['alphafold_program_dir'])

                errormsg = ""

                if not os.path.exists(alndir):
                    errormsg = errormsg + f"Cannot find alignment directory for {targetname}: {alndir}\n"

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

                if len(errormsg) == 0:
                    if not complete_result(method_out_dir, 5 * int(self.params['num_monomer_predictions_per_model'])):
                        temp_struct_csv = os.path.join(template_dir, "sequence_templates.csv")
                        struct_atom_dir = os.path.join(template_dir, "templates")
                        cmd = f"python {self.params['alphafold_program']} " \
                            f"--bfd_uniref_a3m={uniref30_a3m} " \
                            f"--bfd_a3m={bfd_a3m} " \
                            f"--mgnify_sto={mgnify_sto} " \
                            f"--uniref90_sto={uniref90_sto} " \
                            f"--temp_struct_csv={temp_struct_csv} " \
                            f"--struct_atom_dir={struct_atom_dir} " \
                            f"--output_dir={method_out_dir} " + common_parameters
                        cmds += [cmd]
                else:
                    print(errormsg)

            elif run_method ==  "colabfold":
                os.chdir(self.params['alphafold_program_dir'])
                errormsg = ""
                uniref90_sto = os.path.join(alndir, targetname + '_uniref90.sto')
                if not os.path.exists(uniref90_sto):
                    errormsg = errormsg + f"Cannot find uniref90 alignment for {targetname}: {uniref90_sto}\n"

                colabfold_a3m = os.path.join(alndir, targetname + '_colabfold.a3m')
                if not os.path.exists(colabfold_a3m):
                    errormsg = errormsg + f"Cannot find rosettafold alignment for {targetname}: {colabfold_a3m}\n"

                if len(errormsg) == 0:
                    if not complete_result(method_out_dir, 5 * int(self.params['num_monomer_predictions_per_model'])):
                        cmd = f"python {self.params['alphafold_program']} " \
                            f"--custom_msa={colabfold_a3m} " \
                            f"--uniref90_sto={uniref90_sto} " \
                            f"--output_dir={method_out_dir} " + common_parameters
                        cmds += [cmd]
                else:
                    print(errormsg)

            elif run_method ==  "colabfold+seq_template":
                os.chdir(self.params['alphafold_program_dir'])
                errormsg = ""
                uniref90_sto = os.path.join(alndir, targetname + '_uniref90.sto')
                if not os.path.exists(uniref90_sto):
                    errormsg = errormsg + f"Cannot find uniref90 alignment for {targetname}: {uniref90_sto}\n"

                colabfold_a3m = os.path.join(alndir, targetname + '_colabfold.a3m')
                if not os.path.exists(colabfold_a3m):
                    errormsg = errormsg + f"Cannot find rosettafold alignment for {targetname}: {colabfold_a3m}\n"

                if len(errormsg) == 0:
                    if not complete_result(method_out_dir, 5 * int(self.params['num_monomer_predictions_per_model'])):
                        temp_struct_csv = os.path.join(template_dir, "sequence_templates.csv")
                        struct_atom_dir = os.path.join(template_dir, "templates")
                        cmd = f"python {self.params['alphafold_program']} " \
                            f"--custom_msa={colabfold_a3m} " \
                            f"--uniref90_sto={uniref90_sto} " \
                            f"--temp_struct_csv={temp_struct_csv} " \
                            f"--struct_atom_dir={struct_atom_dir} " \
                            f"--output_dir={method_out_dir} " + common_parameters
                        cmds += [cmd]
                else:
                    print(errormsg)
        
            elif run_method ==  "img":
                os.chdir(self.params['alphafold_program_dir'])
                errormsg = ""
                uniref90_sto = os.path.join(alndir, targetname + '_uniref90.sto')
                if not os.path.exists(uniref90_sto):
                    errormsg = errormsg + f"Cannot find uniref90 alignment for {targetname}: {uniref90_sto}\n"

                img_a3m = os.path.join(alndir, targetname + '.a3m')
                if not os.path.exists(img_a3m):
                    errormsg = errormsg + f"Cannot find img alignment for {targetname}: {img_a3m}\n"

                if len(errormsg) == 0:
                    if not complete_result(method_out_dir, 5 * int(self.params['num_monomer_predictions_per_model'])):
                        cmd = f"python {self.params['alphafold_program']} " \
                            f"--custom_msa={img_a3m} " \
                            f"--uniref90_sto={uniref90_sto} " \
                            f"--output_dir={method_out_dir} " + common_parameters 
                        cmds += [cmd]
                else:
                    print(errormsg)

            elif run_method ==  "img+seq_template":
                os.chdir(self.params['alphafold_program_dir'])
                errormsg = ""
                uniref90_sto = os.path.join(alndir, targetname + '_uniref90.sto')
                if not os.path.exists(uniref90_sto):
                    errormsg = errormsg + f"Cannot find uniref90 alignment for {targetname}: {uniref90_sto}\n"

                img_a3m = os.path.join(alndir, targetname + '.a3m')
                if not os.path.exists(img_a3m):
                    errormsg = errormsg + f"Cannot find img alignment for {targetname}: {img_a3m}\n"

                if len(errormsg) == 0:
                    if not complete_result(method_out_dir, 5 * int(self.params['num_monomer_predictions_per_model'])):
                        temp_struct_csv = os.path.join(template_dir, "sequence_templates.csv")
                        struct_atom_dir = os.path.join(template_dir, "templates")
                        cmd = f"python {self.params['alphafold_program']} " \
                            f"--custom_msa={img_a3m} " \
                            f"--uniref90_sto={uniref90_sto} " \
                            f"--temp_struct_csv={temp_struct_csv} " \
                            f"--output_dir={method_out_dir} " + common_parameters
                        cmds += [cmd]
                else:
                    print(errormsg)
            
            elif run_method ==  "dhr":
                os.chdir(self.params['alphafold_program_dir'])
                errormsg = ""
                uniref90_sto = os.path.join(alndir, targetname + '_uniref90.sto')
                if not os.path.exists(uniref90_sto):
                    errormsg = errormsg + f"Cannot find uniref90 alignment for {targetname}: {uniref90_sto}\n"

                dhr_af_a3m = os.path.join(alndir, targetname + '_dhr_af.a3m')
                if not os.path.exists(dhr_af_a3m):
                    errormsg = errormsg + f"Cannot find rosettafold alignment for {targetname}: {dhr_af_a3m}\n"

                if len(errormsg) == 0:
                    if not complete_result(method_out_dir, 5 * int(self.params['num_monomer_predictions_per_model'])):
                        cmd = f"python {self.params['alphafold_program']} " \
                            f"--custom_msa={dhr_af_a3m} " \
                            f"--uniref90_sto={uniref90_sto} " \
                            f"--output_dir={method_out_dir} " + common_parameters
                        cmds += [cmd]
                else:
                    print(errormsg)

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
                os.chdir(self.params['deepfold_program_dir'])
                os.makedirs(method_out_dir, exist_ok=True)
                pickle_path = os.path.join(outdir, self.method2dir['default'], 'features.pkl')
                # if not complete_result(method_out_dir, 5 * int(self.params['num_monomer_predictions_per_model'])):
                if not complete_result(method_out_dir, 5):
                    cmd = f"python {self.params['deepfold_program']} " \
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

            elif run_method.find('deepmsa_noimg') >= 0:

                os.chdir(self.params['alphafold_program_dir'])
                errormsg = ""
                uniref90_sto = os.path.join(alndir, targetname + '_uniref90.sto')
                if not os.path.exists(uniref90_sto):
                    errormsg = errormsg + f"Cannot find uniref90 alignment for {targetname}: {uniref90_sto}\n"

                deepmsa_noimg_a3m = os.path.join(alndir, 'DeepMSA2_a3m', 'MSA', run_method.replace('deepmsa_noimg_', '') + 'a3m')
                if not os.path.exists(deepmsa_noimg_a3m):
                    deepmsa_noimg_a3m = os.path.join(alndir, 'DeepMSA2_a3m', 'MSA', run_method.replace('deepmsa_noimg_', '') + '.a3m')
                    if not os.path.exists(deepmsa_noimg_a3m):
                        errormsg = errormsg + f"Cannot find img alignment for {targetname}: {deepmsa_noimg_a3m}\n"

                if len(errormsg) == 0:
                    if not complete_result(method_out_dir, 5 * int(self.params['num_monomer_predictions_per_model'])):
                        cmd = f"python {self.params['alphafold_program']} " \
                            f"--custom_msa={deepmsa_noimg_a3m} " \
                            f"--uniref90_sto={uniref90_sto} " \
                            f"--output_dir={method_out_dir} " + common_parameters 
                        cmds += [cmd]
                else:
                    print(errormsg)
            
            elif run_method.find('deepmsa_img') >= 0:

                os.chdir(self.params['alphafold_program_dir'])
                errormsg = ""
                uniref90_sto = os.path.join(alndir, targetname + '_uniref90.sto')
                if not os.path.exists(uniref90_sto):
                    errormsg = errormsg + f"Cannot find uniref90 alignment for {targetname}: {uniref90_sto}\n"

                deepmsa_img_a3m = os.path.join(alndir, 'DeepMSA2_a3m', 'JGI', run_method.replace('deepmsa_img_', '') + 'a3m')
                if not os.path.exists(deepmsa_img_a3m):
                    deepmsa_img_a3m = os.path.join(alndir, 'DeepMSA2_a3m', 'JGI', run_method.replace('deepmsa_img_', '') + '.a3m')
                    if not os.path.exists(deepmsa_img_a3m):
                        errormsg = errormsg + f"Cannot find img alignment for {targetname}: {deepmsa_img_a3m}\n"

                if len(errormsg) == 0:
                    if not complete_result(method_out_dir, 5 * int(self.params['num_monomer_predictions_per_model'])):
                        cmd = f"python {self.params['alphafold_program']} " \
                            f"--custom_msa={deepmsa_img_a3m} " \
                            f"--uniref90_sto={uniref90_sto} " \
                            f"--output_dir={method_out_dir} " + common_parameters 
                        cmds += [cmd]
                else:
                    print(errormsg)

            for cmd in cmds:
                try:
                    print(cmd)
                    os.system(cmd)
                except Exception as e:
                    print(e)


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
