# Copyright 2022 Huawei Technologies Co., Ltd
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
# ============================================================================
"""eval script"""
import argparse
import os
import stat
import json
import time
import ast
import numpy as np
import pickle

import mindspore.context as context
import mindspore.common.dtype as mstype
from mindspore import Tensor, nn, save_checkpoint, load_checkpoint, load_param_into_net
from mindsponge.cell.amp import amp_convert
from mindsponge.cell.mask import LayerNormProcess
from mindsponge.common.config_load import load_config
from mindsponge.common.protein import to_pdb, from_prediction

from data import Feature, RawFeatureGenerator, create_dataset, get_crop_size, get_raw_feature, process_pdb
from model import MegaFold, compute_confidence, MegaEvogen
from model.assessment import CombineModel, load_weights
from module.fold_wrapcell import TrainOneStepCell, WithLossCell, WithLossCellAssessment
from module.evogen_block import absolute_position_embedding
from module.lr import cos_decay_lr

parser = argparse.ArgumentParser(description='Inputs for eval.py')
parser.add_argument('--data_config', default="./config/data.yaml", help='data process config')
parser.add_argument('--model_config', default="./config/model.yaml", help='model config')
parser.add_argument('--evogen_config', default="./config/evogen.yaml", help='evogen config')
parser.add_argument('--input_path', help='processed raw feature path')
parser.add_argument('--pdb_path', type=str, help='Location of training pdb file.')
parser.add_argument('--use_pkl', type=ast.literal_eval, default=False,
                    help="use pkl as input or fasta file as input, in default use fasta")
parser.add_argument('--checkpoint_path', help='checkpoint path')
parser.add_argument('--checkpoint_path_assessment', help='assessment model checkpoint path')
parser.add_argument('--device_id', default=0, type=int, help='DEVICE_ID')
parser.add_argument('--is_training', type=ast.literal_eval, default=False, help='is training or not')
parser.add_argument('--run_platform', default='Ascend', type=str, help='which platform to use, Ascend or GPU')
parser.add_argument('--run_distribute', type=ast.literal_eval, default=False, help='run distribute')
parser.add_argument('--resolution_data', type=str, default=None, help='Location of resolution data file.')
parser.add_argument('--loss_scale', type=float, default=1024.0, help='loss scale')
parser.add_argument('--gradient_clip', type=float, default=0.1, help='gradient clip value')
parser.add_argument('--total_steps', type=int, default=9600000, help='total steps')
parser.add_argument('--decoy_pdb_path', type=str, help='Location of decoy pdb file.')
parser.add_argument('--run_assessment', type=int, default=0, help='Run pdb assessment.')
parser.add_argument('--run_evogen', type=int, default=0, help='Run pdb assessment.')
parser.add_argument('--output_dir', type=str,  help='Run pdb assessment.')
arguments = parser.parse_args()


def fold_infer(args):
    '''mega fold inference'''
    data_cfg = load_config(args.data_config)
    model_cfg = load_config(args.model_config)
    data_cfg.eval.crop_size = get_crop_size(args.input_path, args.use_pkl)
    model_cfg.seq_length = data_cfg.eval.crop_size
    slice_key = "seq_" + str(model_cfg.seq_length)
    
    while slice_key not in vars(model_cfg.slice):
        data_cfg.eval.crop_size = data_cfg.eval.crop_size + 256
        model_cfg.seq_length = data_cfg.eval.crop_size
        slice_key = "seq_" + str(model_cfg.seq_length)

    slice_val = vars(model_cfg.slice)[slice_key]
    model_cfg.slice = slice_val

    megafold = MegaFold(model_cfg, mixed_precision=args.mixed_precision)
    load_checkpoint(args.checkpoint_path, megafold)
    if args.mixed_precision:
        fp32_white_list = (nn.Softmax, nn.LayerNorm)
        amp_convert(megafold, fp32_white_list)
    else:
        megafold.to_float(mstype.float32)

    seq_files = []
    
    if not args.use_pkl:
        feature_generator = RawFeatureGenerator(database_search_config=data_cfg.database_search)
        seq_files = [fastafile for fastafile in os.listdir(args.input_path) if fastafile.find('.fasta') > 0]
    else:
        feature_generator = None
        seq_files = [pklfile for pklfile in os.listdir(args.input_path) if pklfile.find('.pkl') > 0]

    for seq_file in seq_files:
        t1 = time.time()
        seq_name = seq_file.split('.')[0]
        os.makedirs(f'{args.output_dir}/{seq_name}', exist_ok=True)
        
        raw_feature = get_raw_feature(os.path.join(args.input_path, seq_file), feature_generator, args.use_pkl)
        # continue
        if not args.use_pkl:
            with open(f'{args.output_dir}/{seq_name}/feature.pkl', 'wb') as f:
                pickle.dump(raw_feature, f, protocol=4)

        ori_res_length = raw_feature['msa'].shape[1]
        processed_feature = Feature(data_cfg, raw_feature)
        feat, prev_pos, prev_msa_first_row, prev_pair = processed_feature.pipeline(data_cfg,
                                                                                   mixed_precision=args.mixed_precision)
        prev_pos = Tensor(prev_pos)
        prev_msa_first_row = Tensor(prev_msa_first_row)
        prev_pair = Tensor(prev_pair)
        t2 = time.time()
        for i in range(data_cfg.common.num_recycle):
            feat_i = [Tensor(x[i]) for x in feat]
            result = megafold(*feat_i,
                              prev_pos,
                              prev_msa_first_row,
                              prev_pair)
            prev_pos, prev_msa_first_row, prev_pair, predicted_lddt_logits = result
        t3 = time.time()
        final_atom_positions = prev_pos.asnumpy()[:ori_res_length]
        final_atom_mask = feat[16][0][:ori_res_length]
        predicted_lddt_logits = predicted_lddt_logits.asnumpy()[:ori_res_length]
        confidence, plddt = compute_confidence(predicted_lddt_logits, return_lddt=True)

        b_factors = plddt[:, None] * final_atom_mask

        unrelaxed_protein = from_prediction(final_atom_positions,
                                            final_atom_mask,
                                            feat[4][0][:ori_res_length],
                                            feat[17][0][:ori_res_length],
                                            b_factors)
        pdb_file = to_pdb(unrelaxed_protein)
        os.makedirs(f'{args.output_dir}/{seq_name}', exist_ok=True)
        os_flags = os.O_RDWR | os.O_CREAT
        os_modes = stat.S_IRWXU
        pdb_path = f'{args.output_dir}/{seq_name}/unrelaxed_{seq_name}.pdb'
        with os.fdopen(os.open(pdb_path, os_flags, os_modes), 'w') as fout:
            fout.write(pdb_file)
        t4 = time.time()
        timings = {"pre_process_time": round(t2 - t1, 2),
                   "predict time ": round(t3 - t2, 2),
                   "pos_process_time": round(t4 - t3, 2),
                   "all_time": round(t4 - t1, 2),
                   "confidence": round(confidence, 2)}

        print(timings)
        with os.fdopen(os.open(f'{args.output_dir}/{seq_name}/timings', os_flags, os_modes), 'w') as fout:
            fout.write(json.dumps(timings))

if __name__ == "__main__":
    if arguments.run_platform == 'Ascend' and not arguments.is_training:
        context.set_context(mode=context.GRAPH_MODE,
                            device_target="Ascend",
                            memory_optimize_level="O1",
                            max_call_depth=6000,
                            device_id=arguments.device_id)
        arguments.mixed_precision = 1
    elif arguments.run_platform == 'Ascend' and arguments.is_training:
        context.set_context(mode=context.GRAPH_MODE,
                            device_target="Ascend",
                            max_device_memory="29GB",
                            device_id=arguments.device_id)
        arguments.mixed_precision = 1
    elif arguments.run_platform == 'GPU':
        context.set_context(mode=context.GRAPH_MODE,
                            device_target="GPU",
                            max_call_depth=6000,
                            graph_kernel_flags="--disable_expand_ops=Softmax --disable_cluster_ops=ReduceSum "
                                               "--composite_op_limit_size=50",
                            device_id=arguments.device_id,
                            enable_graph_kernel=True)
        arguments.mixed_precision = 0
    else:
        raise Exception("Only support GPU or Ascend")

    fold_infer(arguments)
