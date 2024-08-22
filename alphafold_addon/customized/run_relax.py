# Copyright 2021 DeepMind Technologies Limited
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Full AlphaFold protein structure prediction script."""
import enum
import json
import os
import pathlib
import pickle
import random
import shutil
import sys
import time
from typing import Any, Dict, Union

from absl import app
from absl import flags
from absl import logging
from alphafold.common import protein
from alphafold.data_custom import custom_params
from alphafold.common import residue_constants
from alphafold.data_custom import pipeline
from alphafold.data_custom import pipeline_multimer
from alphafold.data_custom import pipeline_multimer_custom
from alphafold.data_custom import templates
from alphafold.data_custom import templates_custom
from alphafold.data_custom.tools import hhsearch
from alphafold.data_custom.tools import hmmsearch
from alphafold.model import config
from alphafold.model import data
from alphafold.model_custom import model
from alphafold.relax import relax
import numpy as np
import jax.numpy as jnp
# Internal import (7716).

USAlign_program = "/bmlfast/bml_casp16/MULTICOM4/tools/USalign"
def cal_tmscore(pdb1, pdb2):
    cmd = USAlign_program + ' ' + pdb1 + ' ' +  pdb2 + " -ter 1 -tmscore 6 | grep TM-score | awk '{print $2}' "
    print(cmd)
    tmscore_contents = os.popen(cmd).read().split('\n')
    tmscore = float(tmscore_contents[1].rstrip('\n'))
    return tmscore

@enum.unique
class ModelsToRelax(enum.Enum):
  ALL = 0
  BEST = 1
  NONE = 2
  TOPN = 3

logging.set_verbosity(logging.INFO)

flags.DEFINE_string('output_dir', None, 'Output directory')
flags.DEFINE_string('database_dir', None, 'Alphafold databse directory')
flags.DEFINE_enum_class('models_to_relax', ModelsToRelax.TOPN, ModelsToRelax,
                        'The models to run the final relaxation step on. '
                        'If `all`, all models are relaxed, which may be time '
                        'consuming. If `best`, only the most confident model '
                        'is relaxed. If `none`, relaxation is not run. Turning '
                        'off relaxation might result in predictions with '
                        'distracting stereochemical violations but might help '
                        'in case you are having issues with the relaxation '
                        'stage.')
flags.DEFINE_boolean('use_gpu_relax', False, 'Whether to relax on GPU. '
                     'Relax on GPU can be much faster than CPU, so it is '
                     'recommended to enable if possible. GPUs must be available'
                     ' if this setting is enabled.')
flags.DEFINE_integer('relax_topn_predictions', 5,
                        'The models to run the final relaxation step on. '
                        'If `all`, all models are relaxed, which may be time '
                        'consuming. If `best`, only the most confident model '
                        'is relaxed. If `none`, relaxation is not run. Turning '
                        'off relaxation might result in predictions with '
                        'distracting stereochemical violations but might help '
                        'in case you are having issues with the relaxation '
                        'stage.')
flags.DEFINE_enum('model_preset', 'multimer',
                  ['monomer', 'monomer_casp14', 'monomer_ptm', 'multimer', 'multimer_v1', 'multimer_v2'],
                  'Choose preset model configuration - the monomer model, '
                  'the monomer model with extra ensembling, monomer model with '
                  'pTM head, or multimer model')

FLAGS = flags.FLAGS

MAX_TEMPLATE_HITS = 20
RELAX_MAX_ITERATIONS = 0
RELAX_ENERGY_TOLERANCE = 2.39
RELAX_STIFFNESS = 10.0
RELAX_EXCLUDE_RESIDUES = []
RELAX_MAX_OUTER_ITERATIONS = 3

def _reorder_chains(pdbstring):
    new_pdbstring = []
    first_chain_id = None
    for line in pdbstring.split('\n'):
        if line.startswith('ATOM') or line.startswith('TER'):
            chain_id = line[21]
            if first_chain_id is None:
                first_chain_id = chain_id
            if first_chain_id != "A":
                new_pdbstring += [line[:21] + protein.PDB_CHAIN_IDS[protein.PDB_CHAIN_IDS.find(chain_id)-1] + line[22:]]
            else:
                new_pdbstring += [line]
        else:
            new_pdbstring += [line]
    return '\n'.join(new_pdbstring)

def main(argv):
    if len(argv) > 1:
        raise app.UsageError('Too many command-line arguments.')

    feature_pickle = os.path.join(FLAGS.output_dir, 'features.pkl')
    with open(feature_pickle, 'rb') as f:
        feature_dict = pickle.load(f)
    # print(feature_dict)

    # amber_relaxer = None
    amber_relaxer = relax.AmberRelaxation(
        max_iterations=RELAX_MAX_ITERATIONS,
        tolerance=RELAX_ENERGY_TOLERANCE,
        stiffness=RELAX_STIFFNESS,
        exclude_residues=RELAX_EXCLUDE_RESIDUES,
        max_outer_iterations=RELAX_MAX_OUTER_ITERATIONS,
        use_gpu=FLAGS.use_gpu_relax)

    model_runners = {}
    model_names = config.MODEL_PRESETS[FLAGS.model_preset]
    for model_name in model_names:
        model_config = config.model_config(model_name)
        model_params = data.get_model_haiku_params(model_name=model_name, data_dir=FLAGS.database_dir)
        model_runner = model.RunModel(model_config, model_params)
        model_runners[model_name] = model_runner

    ranking_confidences = {}
    ranked_order = []
    for pdb in os.listdir(FLAGS.output_dir):
        if pdb.find('unrelaxed') != 0:
            continue
        pkl = pdb.replace('unrelaxed', 'result').replace('.pdb', '.pkl')
        with open(FLAGS.output_dir + '/' + pkl, 'rb') as f:
            prediction_result = pickle.load(f)
            print(prediction_result['ranking_confidence'])
            ranking_confidences[pdb.replace('unrelaxed_', '').replace('.pdb', '')] = float(prediction_result['ranking_confidence'])
            ranked_order.append(pdb.replace('unrelaxed_', '').replace('.pdb', ''))

    # Rank by model confidence and write out relaxed PDBs in rank order.
    ranked_order = []
    for idx, (model_name, _) in enumerate(
            sorted(ranking_confidences.items(), key=lambda x: x[1], reverse=True)):
        ranked_order.append(model_name)
        ranked_output_path = os.path.join(FLAGS.output_dir, f'ranked_{idx}.pdb')
        srcpdb = os.path.join(FLAGS.output_dir, f"relaxed_{model_name}.pdb")
        if not os.path.exists(srcpdb):
            if idx >= FLAGS.relax_topn_predictions:
                srcpdb = os.path.join(FLAGS.output_dir, f"unrelaxed_{model_name}.pdb")
                pdbstring = ''.join(open(srcpdb).readlines())
                with open(ranked_output_path, 'w') as f:
                    f.write(_reorder_chains(pdbstring))
            else:
                random_seed = random.randrange(sys.maxsize // len(model_runners))
                logging.info('Using random seed %d for the data pipeline', random_seed)
                for name in model_runners:
                    if model_name.find(name) == 0:
                        processed_feature_dict = model_runners[name].process_features(feature_dict, random_seed=random_seed)
                        break

                result_pickle = os.path.join(FLAGS.output_dir, f"result_{model_name}.pkl")
                with open(result_pickle, 'rb') as f:
                    prediction_result=pickle.load(f)
                    plddt = prediction_result['plddt']
                    ranking_confidences[model_name] = prediction_result['ranking_confidence']
                    # Add the predicted LDDT in the b-factor column.
                    # Note that higher predicted LDDT value means higher model confidence.
                    plddt_b_factors = np.repeat(
                        plddt[:, None], residue_constants.atom_type_num, axis=-1)
                    unrelaxed_protein = protein.from_prediction(
                        features=processed_feature_dict,
                        result=prediction_result,
                        b_factors=plddt_b_factors,
                        remove_leading_feature_dimension=not model_runner.multimer_mode)

                    t_0 = time.time()
                    relaxed_pdb_str, _, _ = amber_relaxer.process(prot=unrelaxed_protein)
                    logging.info(f'Relax took {time.time()-t_0} s')
                    # Save the relaxed PDB.
                    with open(srcpdb, 'w') as f:
                        f.write(relaxed_pdb_str)
                    os.system(f"cp {srcpdb} {ranked_output_path}")

        else:
            os.system(f"cp {srcpdb} {ranked_output_path}")

        if idx < FLAGS.relax_topn_predictions:
            unrelaxed_pdb = os.path.join(FLAGS.output_dir, f"unrelaxed_{model_name}.pdb")
            relaxed_pdb = os.path.join(FLAGS.output_dir, f"relaxed_{model_name}.pdb")
            tmscore = cal_tmscore(unrelaxed_pdb, relaxed_pdb)
            if tmscore < 0.9:
                print(f"Warning!!!!!!!!!!! TMscore < 0.9 ({tmscore}) between relaxed and unrelaxed model!")
            else:
                 print(f"TMscore >= 0.9 ({tmscore}) between relaxed and unrelaxed model!")

    ranking_output_path = os.path.join(FLAGS.output_dir, 'ranking_debug.json')
    with open(ranking_output_path, 'w') as f:
        label = 'iptm+ptm' if 'iptm' in prediction_result else 'plddts'
        f.write(json.dumps(
            {label: ranking_confidences, 'order': ranked_order}, indent=4))

if __name__ == '__main__':
    flags.mark_flags_as_required([
        'database_dir',
        'output_dir'
    ])
    app.run(main)
