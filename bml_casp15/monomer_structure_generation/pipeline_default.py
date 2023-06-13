import copy
import os
import sys
import time
from bml_casp15.common.util import makedir_if_not_exists, check_dirs
import pandas as pd
from multiprocessing import Pool
import pathlib
from bml_casp15.common.protein import complete_result


class Monomer_structure_prediction_pipeline_default:

    def __init__(self, params):

        self.params = params

    def process_single(self, fasta_path, alndir, outdir):

        targetname = pathlib.Path(fasta_path).stem

        cmd = ""

        makedir_if_not_exists(outdir)

        os.chdir(self.params['alphafold_program_dir'])

        errormsg = ""

        if not os.path.exists(alndir):
            errormsg = errormsg + f"Cannot find alignment directory for {targetname}: {alndir}\n"

        bfd_uniref30_a3m = alndir + '/' + targetname + '_uniref30_bfd.a3m'
        if not os.path.exists(bfd_uniref30_a3m):
            errormsg = errormsg + f"Cannot find uniclust30 alignment for {targetname}: {bfd_uniref30_a3m}\n"

        mgnify_sto = alndir + '/' + targetname + '_mgnify.sto'
        if not os.path.exists(mgnify_sto):
            errormsg = errormsg + f"Cannot find mgnify alignment for {targetname}: {mgnify_sto}\n"

        uniref90_sto = alndir + '/' + targetname + '_uniref90.sto'
        if not os.path.exists(uniref90_sto):
            errormsg = errormsg + f"Cannot find uniref90 alignment for {targetname}: {uniref90_sto}\n"

        if len(errormsg) == 0:
            if not complete_result(f"{outdir}/default", 5 * int(self.params['num_monomer_predictions_per_model'])):
                try:
                    cmd = f"python {self.params['alphafold_default_program']} " \
                          f"--fasta_path {fasta_path} " \
                          f"--env_dir {self.params['alphafold_env_dir']} " \
                          f"--database_dir {self.params['alphafold_database_dir']} " \
                          f"--bfd_uniref_a3ms {bfd_uniref30_a3m} " \
                          f"--mgnify_stos {mgnify_sto} " \
                          f"--uniref90_stos {uniref90_sto} " \
                          f"--monomer_num_ensemble {self.params['monomer_num_ensemble']} " \
                          f"--monomer_num_recycle {self.params['monomer_num_recycle']} " \
                          f"--num_monomer_predictions_per_model {self.params['num_monomer_predictions_per_model']} " \
                          f"--models_to_relax=best " \
                          f"--output_dir {outdir}/default"
                    print(cmd)
                    os.system(cmd)
                except Exception as e:
                    print(e)
        else:
            print(errormsg)

    def process(self, monomers, alndir, outdir, templatedir=None):
        outdir = os.path.abspath(outdir) + "/"
        for fasta_path in monomers:
            fasta_name = pathlib.Path(fasta_path).stem
            monomer_aln_dir = alndir + '/' + fasta_name
            monomer_outdir = outdir + '/' + fasta_name
            monomer_template_dir = ""
            if templatedir is not None:
                monomer_template_dir = templatedir + '/' + fasta_name
            self.process_single(fasta_path=fasta_path, alndir=monomer_aln_dir, outdir=monomer_outdir,
                                template_dir=monomer_template_dir)

        print("The tertiary structure generation for monomers has finished!")
