import copy
import os
import sys
import time, json
from multicom4.common.util import makedir_if_not_exists, check_dirs
import dataclasses
from multicom4.monomer_structure_refinement.iterative_refine_pipeline_v4_50 import *
import pandas as pd
from multicom4.monomer_structure_refinement.util import cal_tmscore
from multicom4.common import config

class refinement_input:
    def __init__(self, fasta_path, pdb_path, pkl_path, msa_path):
        self.fasta_path = fasta_path
        self.pdb_path = pdb_path
        self.pkl_path = pkl_path
        self.msa_path = msa_path


class Monomer_iterative_refinement_pipeline_server:

    def __init__(self, params, config_name):
        self.params = params
        self.config_name = config_name

    def search(self, refinement_inputs, outdir, uniref90_sto=""):
        result_dirs = []

        pipeline = Monomer_iterative_refinement_pipeline(self.params, self.config_name)

        for refine_param in refinement_inputs:
            result_dir = pipeline.search_single(fasta_path=refine_param.fasta_path, pdb_path=refine_param.pdb_path,
                                                pkl_path=refine_param.pkl_path, msa_path=refine_param.msa_path,
                                                outdir=os.path.join(outdir, pathlib.Path(refine_param.pdb_path).stem),
                                                uniref90_sto=uniref90_sto)
            result_dirs += [result_dir]

        return result_dirs


class Monomer_refinement_model_selection(config.pipeline):

    def __init__(self, params, config_name):

        super().__init__()

        self.params = params

        self.predictor_config = self.monomer_config.predictors.foldseek_refine

    def select_v1(self, indir, outdir):
        if os.path.exists(outdir):
            os.system(f"rm -rf {outdir}")
        makedir_if_not_exists(outdir)

        for pdb in os.listdir(indir):
            iteration_dir = os.path.join(indir, pdb, 'iteration1')
            if not os.path.exists(iteration_dir):
                continue
                
            # start_pdb = os.path.join(iteration_dir, 'start.pdb')
            # start_pkl = os.path.join(iteration_dir, 'start.pkl') 
            # start_a3m = os.path.join(iteration_dir, 'start.a3m') 
            # os.system("cp " + start_pdb + " " + os.path.join(outdir, pdb + "_ori.pdb"))
            # os.system("cp " + start_pkl + " " + os.path.join(outdir, pdb + "_ori.pkl"))
            # os.system("cp " + start_a3m + " " + os.path.join(outdir, pdb + "_ori.a3m"))

            refine_pdb = os.path.join(indir, pdb, 'final', 'final.pdb')
            refine_pkl = os.path.join(indir, pdb, 'final', 'final.pkl')
            refine_a3m = os.path.join(indir, pdb, 'final', 'final.a3m')
            os.system("cp " + refine_pdb + " " + os.path.join(outdir, pdb + "_ref.pdb"))
            os.system("cp " + refine_pkl + " " + os.path.join(outdir, pdb + "_ref.pkl"))
            os.system("cp " + refine_a3m + " " + os.path.join(outdir, pdb + "_ref.a3m"))

        pdbs = []
        plddts = []
        for pkl in os.listdir(outdir):
            if pkl.find('.pkl') > 0:
                with open(os.path.join(outdir, pkl), 'rb') as f:
                    prediction_result = pickle.load(f)
                    pdbs += [pkl.replace('.pkl', '.pdb')]
                    plddts += [np.mean(prediction_result['plddt'])]

        df = pd.DataFrame({'model': pdbs, 'plddt': plddts})
        df = df.sort_values(by=['plddt'], ascending=False)
        df.reset_index(inplace=True, drop=True)
        df.to_csv(os.path.join(outdir, 'final_ranking.csv'))

        # if prefix == "refine":
        #     selected_models = []
        #     for i in range(4):
        #         pdb_name = df.loc[i, 'model']
        #         prefix_name = f"{prefix}{i+1}"
        #         os.system("cp " + os.path.join(outdir, pdb_name) + " " + os.path.join(outdir, prefix_name+".pdb"))
        #         os.system("cp " + os.path.join(outdir, pdb_name.replace('.pdb', '.pkl')) + " " + os.path.join(outdir, prefix_name+".pkl"))
        #         os.system("cp " + os.path.join(outdir, pdb_name.replace('.pdb', '.a3m')) + " " + os.path.join(outdir, prefix_name+".a3m"))
        #         selected_models += [pdb_name]

        #     top1_model = df.loc[0, 'model']
        #     added = False
        #     for i in range(4, len(df)):
        #         pdb_name = df.loc[i, 'model']
        #         tmscore, gdtscore = cal_tmscore(self.params['tmscore_program'],
        #                                         os.path.join(outdir, pdb_name),
        #                                         os.path.join(outdir, top1_model),
        #                                         os.path.join(outdir, 'tmp'))
        #         if tmscore < 0.98:
        #             os.system("cp " + os.path.join(outdir, pdb_name) + " " + os.path.join(outdir, prefix+"5.pdb"))
        #             os.system("cp " + os.path.join(outdir, pdb_name.replace('.pdb', '.pkl')) + " " + os.path.join(outdir, prefix+"5.pkl"))
        #             os.system("cp " + os.path.join(outdir, pdb_name.replace('.pdb', '.a3m')) + " " + os.path.join(outdir, prefix+"5.a3m"))
        #             added = True
        #             selected_models += [pdb_name]
        #             break
        #         else:
        #             print(f"The tmscore between {pdb_name} and {top1_model} is larger than 0.98 ({tmscore}), skipped!")
        #     if not added:
        #         pdb_name = df.loc[4, 'model']
        #         os.system("cp " + os.path.join(outdir, pdb_name) + " " + os.path.join(outdir, prefix+"5.pdb"))
        #         os.system("cp " + os.path.join(outdir, pdb_name.replace('.pdb', '.pkl')) + " " + os.path.join(outdir, prefix+"5.pkl"))
        #         os.system("cp " + os.path.join(outdir, pdb_name.replace('.pdb', '.a3m')) + " " + os.path.join(outdir, prefix+"5.a3m"))
        #         selected_models += [pdb_name]

        #     selected_df = pd.DataFrame({'selected_models': selected_models})
        #     selected_df.to_csv(os.path.join(outdir, f'{prefix}_selected.csv'))

        # else:
        #     selected_models = []
        #     for i in range(5):
        #         pdb_name = df.loc[i, 'model']
        #         prefix_name = f"{prefix}{i+1}"
        #         os.system("cp " + os.path.join(outdir, pdb_name) + " " + os.path.join(outdir, prefix_name+".pdb"))
        #         os.system("cp " + os.path.join(outdir, pdb_name.replace('.pdb', '.pkl')) + " " + os.path.join(outdir, prefix_name+".pkl"))
        #         os.system("cp " + os.path.join(outdir, pdb_name.replace('.pdb', '.a3m')) + " " + os.path.join(outdir, prefix_name+".a3m"))
        #         selected_models += [pdb_name]
        #     selected_df = pd.DataFrame({'selected_models': selected_models})
        #     selected_df.to_csv(os.path.join(outdir, f'{prefix}_selected.csv'))

        # selected_models = []
        # for i in range(self.predictor_config.number_of_output_models):
        #     pdb_name = df.loc[i, 'model']
        #     prefix_name = f"{prefix}{i+1}"
        #     os.system("cp " + os.path.join(outdir, pdb_name) + " " + os.path.join(outdir, prefix_name+".pdb"))
        #     os.system("cp " + os.path.join(outdir, pdb_name.replace('.pdb', '.pkl')) + " " + os.path.join(outdir, prefix_name+".pkl"))
        #     os.system("cp " + os.path.join(outdir, pdb_name.replace('.pdb', '.a3m')) + " " + os.path.join(outdir, prefix_name+".a3m"))
        #     selected_models += [pdb_name]
        # selected_df = pd.DataFrame({'selected_models': selected_models})
        # selected_df.to_csv(os.path.join(outdir, f'{prefix}_selected.csv'))

        return outdir

    # def select_v2(self, ranking_df, indir, outdir, prefix, refpdb):
    #     if os.path.exists(outdir):
    #         os.system(f"rm -rf {outdir}")
    #     makedir_if_not_exists(outdir)

    #     start_pdbs = []
    #     refine_pdbs = []
    #     start_plddts = []
    #     refine_plddts = []
    #     final_models = []
    #     tmscores = []
    #     reftmscores = []
    #     for i in range(5):
    #         pdb_name = ranking_df.loc[i, 'model'].replace('.pdb', '')
    #         iteration_dir = os.path.join(indir, pdb_name, 'iteration1')
    #         if not os.path.exists(iteration_dir):
    #             continue
    #         start_pdb = os.path.join(iteration_dir, 'start.pdb')
    #         start_pkl = os.path.join(iteration_dir, 'start.pkl') 
    #         os.system("cp " + start_pdb + " " + os.path.join(outdir, pdb_name + "_ori.pdb"))
    #         os.system("cp " + start_pkl + " " + os.path.join(outdir, pdb_name + "_ori.pkl"))
    #         with open(start_pkl, 'rb') as f:
    #             plddt_start = np.mean(pickle.load(f)['plddt'])

    #         start_pdbs += [f"{pdb_name}_ori.pdb"]
    #         start_plddts += [plddt_start]

    #         refine_pdb = os.path.join(indir, pdb_name, 'final', 'final.pdb')
    #         refine_pkl = os.path.join(indir, pdb_name, 'final','final.pkl')
    #         os.system("cp " + refine_pdb + " " + os.path.join(outdir, pdb_name + "_ref.pdb"))
    #         os.system("cp " + refine_pkl + " " + os.path.join(outdir, pdb_name + "_ref.pkl"))
    #         with open(refine_pkl, 'rb') as f:
    #             plddt_ref = np.mean(pickle.load(f)['plddt'])

    #         refine_pdbs += [f"{pdb_name}_ref.pdb"]
    #         refine_plddts += [plddt_ref]
    #         tmscore, _ = cal_tmscore(tmscore_program=self.params['tmscore_program'],
    #                                  inpdb=start_pdb,
    #                                  nativepdb=refine_pdb,
    #                                  tmpdir=os.path.join(outdir, 'tmp'))
    #         tmscores += [tmscore]

    #         reftmscore = 0
    #         if os.path.exists(refpdb):
    #             reftmscore, _ = cal_tmscore(tmscore_program=self.params['tmscore_program'],
    #                                      inpdb=refpdb,
    #                                      nativepdb=refine_pdb,
    #                                      tmpdir=os.path.join(outdir, 'tmp'))
    #         reftmscores += [reftmscore]

    #         prefix_name = f"{prefix}{i+1}"
                
    #         if plddt_start > plddt_ref:
    #             os.system("cp " + os.path.join(outdir, pdb_name + "_ori.pdb") + " " + os.path.join(outdir, prefix_name+".pdb"))
    #             os.system("cp " + os.path.join(outdir, pdb_name + "_ori.pkl") + " " + os.path.join(outdir, prefix_name+".pkl"))
    #             final_models += [f"{pdb_name}_ori.pdb"]
    #         else:
    #             os.system("cp " + os.path.join(outdir, pdb_name + "_ref.pdb") + " " + os.path.join(outdir, prefix_name+".pdb"))
    #             os.system("cp " + os.path.join(outdir, pdb_name + "_ref.pkl") + " " + os.path.join(outdir, prefix_name+".pkl"))
    #             final_models += [f"{pdb_name}_ref.pdb"]

    #     df = pd.DataFrame({'start_model': start_pdbs, 'start_plddt': start_plddts,
    #                        'refine_model': refine_pdbs, 'refine_plddt': refine_plddts,
    #                        'tmscore': tmscores, 'reftmscore': reftmscores, 'final_model': final_models})
    #     df.to_csv(os.path.join(outdir, 'final_ranking.csv'))

    #     return outdir

    def make_predictor_results(self, indir, outdir):

        ranking_file = indir + '/final_ranking.csv'
        ranking_df = pd.read_csv(ranking_file)
        ranking_confidences = {}
        ranked_order = []
        model_count = 0
        for i in range(len(ranking_df)):
            model_name = ranking_df.loc[i, 'model']
            if not model_name.endswith('_ref.pdb'):
                continue

            if not os.path.exists(f"{indir}/{model_name}"):
                raise Exception(f"Cannot find {indir}/{model_name}")

            trg_model_name = f"ranked_{model_count}"

            result_output_path = f"{indir}/{model_name.replace('.pdb', '.pkl')}"
            with open(result_output_path, 'rb') as f:
                prediction_result = pickle.load(f)
                ranking_confidences[trg_model_name] = prediction_result['ranking_confidence']
                ranked_order.append(trg_model_name)
            os.system(f"cp {indir}/{model_name} {outdir}/{trg_model_name}.pdb")
            os.system(f"cp {result_output_path} {outdir}/result_{trg_model_name}.pkl")
            model_count += 1

        ranking_output_path = os.path.join(outdir, 'ranking_debug.json')
        with open(ranking_output_path, 'w') as f:
            label = 'iptm+ptm' if 'iptm' in prediction_result else 'plddts'
            f.write(json.dumps(
                {label: ranking_confidences, 'order': ranked_order}, indent=4))

        msadir = outdir + '/msas'
        makedir_if_not_exists(msadir)
        os.system(f"cp {indir}/{model_name.replace('.pdb', '.pkl')} {msadir}/monomer_final.a3m")


