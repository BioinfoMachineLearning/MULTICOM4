import copy
import os
import sys
import time
from multicom4.common.util import makedir_if_not_exists, check_dirs
from multicom4.common.protein import complete_result, parse_fasta
import pandas as pd
from multiprocessing import Pool
import pathlib


class Multimer_structure_prediction_pipeline_default:

    def __init__(self, params):  # , run_methods):

        self.params = params

    def process(self,
                fasta_path,
                chain_id_map,
                aln_dir,
                output_dir):

        makedir_if_not_exists(output_dir)

        common_parameters = f"--fasta_path={fasta_path} " \
                            f"--env_dir={self.params['alphafold_env_dir']} " \
                            f"--database_dir={self.params['alphafold_database_dir']} " \
                            f"--multimer_num_ensemble={self.params['multimer_num_ensemble']} " \
                            f"--multimer_num_recycle={self.params['multimer_num_recycle']} " \
                            f"--num_multimer_predictions_per_model={self.params['num_multimer_predictions_per_model']} " \
                            f"--model_preset={self.params['multimer_model_preset']} " \
                            f"--benchmark={self.params['alphafold_benchmark']} " \
                            f"--use_gpu_relax={self.params['use_gpu_relax']} " \
                            f"--models_to_relax={self.params['models_to_relax']} " \
                            f"--max_template_date={self.params['max_template_date']} "   

        # run alphafold default pipeline:
        outdir = f"{output_dir}/default_multimer"
        monomers = [chain_id for chain_id in chain_id_map]
        if not complete_result(outdir, 5 * int(self.params['num_multimer_predictions_per_model'])):
            os.chdir(self.params['alphafold_program_dir'])
            bfd_uniclust_a3ms = []
            mgnify_stos = []
            uniref90_stos = []
            uniprot_stos = []
            for chain_id in chain_id_map:
                monomer = chain_id
                monomer_bfd_uniclust_a3m = f"{aln_dir}/{monomer}/{monomer}_uniref30_bfd.a3m"
                if not os.path.exists(monomer_bfd_uniclust_a3m):
                    raise Exception(f"Cannot find bfd and uniclust a3m for {monomer}: {monomer_bfd_uniclust_a3m}")
                bfd_uniclust_a3ms += [monomer_bfd_uniclust_a3m]
        
                monomer_mgnify_sto = f"{aln_dir}/{monomer}/{monomer}_mgnify.sto"
                if not os.path.exists(monomer_mgnify_sto):
                    raise Exception(f"Cannot find mgnify sto for {monomer}: {monomer_mgnify_sto}")
                mgnify_stos += [monomer_mgnify_sto]
        
                monomer_uniref90_sto = f"{aln_dir}/{monomer}/{monomer}_uniref90.sto"
                if not os.path.exists(monomer_uniref90_sto):
                    raise Exception(f"Cannot find uniref90 sto for {monomer}: {monomer_uniref90_sto}")
                uniref90_stos += [monomer_uniref90_sto]
        
                monomer_uniprot_sto = f"{aln_dir}/{monomer}/{monomer}_uniprot.sto"
                if not os.path.exists(monomer_uniprot_sto):
                    raise Exception(f"Cannot find uniprot sto for {monomer}: {monomer_uniprot_sto}")
                uniprot_stos += [monomer_uniprot_sto]
        
            cmd = f"python {self.params['alphafold_default_program']} " \
                  f"--bfd_uniref_a3ms {','.join(bfd_uniclust_a3ms)} " \
                  f"--mgnify_stos {','.join(mgnify_stos)} " \
                  f"--uniref90_stos {','.join(uniref90_stos)} " \
                  f"--uniprot_stos {','.join(uniprot_stos)} " \
                  f"--output_dir {outdir} " + common_parameters
        
            print(cmd)
            os.system(cmd)
