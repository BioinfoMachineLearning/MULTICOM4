import os, sys, argparse, time
from multiprocessing import Pool
from multicom_dev.common.util import check_file, check_dir, check_dirs, makedir_if_not_exists, check_contents, \
    read_option_file
from multicom_dev.monomer_alignment_generation.alignment import read_fasta, write_fasta
from multicom_dev.monomer_alignment_generation.pipeline import Monomer_alignment_generation_pipeline
from multicom_dev.monomer_alignments_concatenation.pipeline import *
from multicom_dev.monomer_structure_generation.pipeline import *
from multicom_dev.monomer_templates_concatenation import sequence_based_pipeline_complex_pdb, \
    sequence_based_pipeline_pdb, sequence_based_pipeline, structure_based_pipeline_v2
from multicom_dev.multimer_structure_generation.pipeline import *
from multicom_dev.multimer_structure_generation.iterative_search_pipeline_v0_2 import *
from multicom_dev.multimer_structure_refinement import iterative_refine_pipeline
from absl import flags
from absl import app
import copy

flags.DEFINE_string('option_file', None, 'option file')
flags.DEFINE_string('fasta_path', None, 'Path to multimer fastas')
flags.DEFINE_string('model_dir', None, 'Output directory')
flags.DEFINE_string('output_dir', None, 'Output directory')
FLAGS = flags.FLAGS

PDB_CHAIN_IDS = 'BCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789'
PDB_MAX_CHAINS = len(PDB_CHAIN_IDS)  # := 62.


class FastaChain:
    sequence: str
    description: str


def _make_chain_id_map(sequences, descriptions):
    """Makes a mapping from PDB-format chain ID to sequence and description."""
    if len(sequences) != len(descriptions):
        raise ValueError('sequences and descriptions must have equal length. '
                         f'Got {len(sequences)} != {len(descriptions)}.')
    if len(sequences) > PDB_MAX_CHAINS:
        raise ValueError('Cannot process more chains than the PDB format supports. '
                         f'Got {len(sequences)} chains.')
    chain_id_map = {}
    chain_id_seq_map = {}
    for chain_id, sequence, description in zip(
            PDB_CHAIN_IDS, sequences, descriptions):
        chain_id_map[chain_id] = FastaChain(sequence=sequence, description=description)
        chain_id_seq_map[chain_id] = sequence
    return chain_id_map, chain_id_seq_map

def read_qa_txt(infile):
    models = []
    scores = []
    for line in open(infile):
        line = line.rstrip('\n')
        contents = line.split()
        if contents[0] == "PFRMAT" or contents[0] == "TARGET" or contents[0] == "MODEL" or contents[0] == "QMODE" or \
                contents[0] == "END":
            continue
        model, score = line.split()
        models += [model]
        scores += [float(score)]
    df = pd.DataFrame({'model': models, 'score': scores})
    df = df.sort_values(by=['score'], ascending=False)
    df.reset_index(inplace=True)
    return df


def parse_fasta(fasta_string):
    sequences = []
    descriptions = []
    index = -1
    for line in fasta_string.splitlines():
        line = line.strip()
        if line.startswith('>'):
            index += 1
            descriptions.append(line[1:])  # Remove the '>' at the beginning.
            sequences.append('')
            continue
        elif not line:
            continue  # Skip blank lines.
        sequences[index] += line
    return sequences, descriptions


def main(argv):
    if len(argv) > 1:
        raise app.UsageError('Too many command-line arguments.')

    check_file(FLAGS.option_file)

    params = read_option_file(FLAGS.option_file)

    makedir_if_not_exists(FLAGS.output_dir)

    check_dirs(params, ['hhblits_program', 'jackhmmer_program'], isdir=False)

    check_file(fasta_path)

    makedir_if_not_exists(FLAGS.output_dir)

    with open(fasta_path) as f:
        input_fasta_str = f.read()
    input_seqs, input_descs = parsers.parse_fasta(input_fasta_str)
    chain_id_map, chain_id_seq_map = _make_chain_id_map(sequences=input_seqs,
                                                        descriptions=input_descs)

    print("8. Start to evaluate multimer models")

    N7_outdir = FLAGS.model_dir

    N8_outdir = outdir + '/N8_multimer_structure_evaluation'
    pipeline = Multimer_structure_evaluation_pipeline(params=params)
    multimer_qa_result = pipeline.process(N7_outdir, N8_outdir)

    print("#################################################################################################")

    print("#################################################################################################")

    # print("9. Start to refine multimer models based on the qa rankings")

    # N9_outdir = outdir + '/N9_multimer_structure_refinement'

    # makedir_if_not_exists(N9_outdir)

    # ref_ranking = read_qa_txt(multimer_qa_result['alphafold'])  # apollo or average ranking or the three qas

    # pipeline = iterative_refine_pipeline.Multimer_iterative_refinement_pipeline_server(params=params)
    # refine_inputs = []
    # for i in range(5):
    #     pdb_name = ref_ranking.loc[i, 'model']
    #     msa_paths = {}
    #     for chain_id in chain_id_map:
    #         msa_paths[chain_id] = dict(
    #             paired_msa=f"{N8_outdir}/msa/{chain_id_map[chain_id].description}/{pdb_name}.multimer.a3m",
    #             monomer_msa=f"{N8_outdir}/msa/{chain_id_map[chain_id].description}/{pdb_name}.monomer.a3m")

    #     refine_input = iterative_refine_pipeline.refinement_input(chain_id_map=chain_id_map,
    #                                                               fasta_path=fasta_path,
    #                                                               pdb_path=N8_outdir + '/pdb/' + pdb_name,
    #                                                               pkl_path=N8_outdir + '/pkl/' + pdb_name.replace(
    #                                                                   '.pdb', '.pkl'),
    #                                                               msa_paths=msa_paths)
    #     refine_inputs += [refine_input]
    # pipeline.search(refinement_inputs=refine_inputs, outdir=N9_outdir)

    # final_dir = N9_outdir + '/final'
    # makedir_if_not_exists(final_dir)

    # pipeline = iterative_refine_pipeline.Multimer_refinement_model_selection()
    # pipeline.select_v1(indir=N9_outdir, outdir=final_dir)

    # print("The refinement for the top-ranked multimer models has been finished!")


if __name__ == '__main__':
    flags.mark_flags_as_required([
        'option_file',
        'fasta_path',
        'output_dir'
    ])
    app.run(main)
