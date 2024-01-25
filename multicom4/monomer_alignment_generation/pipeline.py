import os, sys, argparse, time
from multiprocessing import Pool
from tqdm import tqdm
from multicom4.common.util import is_file, is_dir, makedir_if_not_exists, check_contents, read_option_file, check_dirs
from multicom4.monomer_alignment_generation.alignment import *
from multicom4.monomer_alignment_generation.rosettafold_msa_runner import *
from multicom4.monomer_alignment_generation.colabfold_msa_runner import *
from multicom4.monomer_alignment_generation.img_msa_runner import *
from multicom4.monomer_alignment_generation.deepmsa2_runner import *
from multicom4.monomer_alignment_generation.dhr_runner import *
from multicom4.tool import hhblits
from multicom4.tool import jackhmmer
import pathlib
from multicom4.monomer_templates_concatenation import parsers

def run_msa_tool(inparams):
    msa_runner, input_fasta_path, msa_out_path, msa_out_name, msa_key = inparams
    """Runs an MSA tool, checking if output already exists first."""
    msa_out_file = os.path.join(msa_out_path, msa_out_name)
    if not os.path.exists(msa_out_file) or len(open(msa_out_file).readlines()) == 0:
        workdir = os.path.join(msa_out_path, msa_key)
        makedir_if_not_exists(workdir)
        if msa_key == "DeepMSA2_a3m":
            result = msa_runner.query(input_fasta_path, workdir)
            os.system(f"cp {result} {msa_out_file}")
        else:
            temp_msa_out_file = os.path.join(workdir, msa_out_name)
            result = msa_runner.query(input_fasta_path, temp_msa_out_file)
            os.system(f"cp {temp_msa_out_file} {msa_out_file}")

    return msa_key, msa_out_file

def combine_a3ms(msas, msa_save_path, max_seq = 100000):
    seen_desc = []
    seen_sequences = []
    for msa_index, msa in enumerate(msas):
        if not msa:
            raise ValueError(f'MSA {msa_index} must contain at least one sequence.')
        for sequence_index, sequence in enumerate(msa.sequences):
            if sequence in seen_sequences:
                continue
            seen_sequences += [sequence]
            seen_desc += [msa.descriptions[sequence_index]]
            if len(seen_sequences) > max_seq:
                break

        if len(seen_sequences) > max_seq:
            break

    with open(msa_save_path, 'w') as fw:
        for (desc, seq) in zip(seen_desc, seen_sequences):
            fw.write(f'>{desc}\n{seq}\n')

class Monomer_alignment_generation_pipeline:
    """Runs the alignment tools and assembles the input features."""

    def __init__(self,
                 jackhmmer_binary_path,
                 hhblits_binary_path,
                 colabfold_search_binary,
                 colabfold_split_msas_binary,
                 mmseq_binary,
                 deepmsa2_path,
                 dhr_program_path,
                 uniref90_database_path,
                 mgnify_database_path,
                 small_bfd_database_path,
                 bfd_database_path,
                 uniref30_database_path,
                 uniclust30_database_path,
                 uniprot_database_path,
                 JGIclust_database_path,
                 metaclust_database_path,
                 colabfold_databases,
                 dhr_database_path,
                 hhfilter_binary_path="",
                 mgnify_max_hits: int = 501,
                 uniref_max_hits: int = 10000,
                 use_precomputed_msas: bool = False):
        """Initializes the data pipeline."""

        # alignment generation pipeline from alphafold
        self.mgnify_max_hits = mgnify_max_hits
        self.uniref_max_hits = uniref_max_hits
        self.hhfilter_binary_path = hhfilter_binary_path
        
        self.jackhmmer_uniref90_runner = None
        self.hhblits_bfd_runner = None
        self.hhblits_uniref_runner = None
        self.jackhmmer_mgnify_runner = None
        self.hhblits_uniclust_runner = None
        self.jackhmmer_uniprot_runner = None
        self.hhblits_uniclust_folddock_runner = None
        self.jackhmmer_small_bfd_runner = None
        self.rosettafold_msa_runner = None
        self.colabfold_msa_runner = None
        self.uniref30_bfd_msa_runner = None
        self.unclust30_bfd_msa_runner = None
        self.deepmsa2_runner = None
        self.dhr_runner = None

        if len(uniref90_database_path) > 0:
            self.jackhmmer_uniref90_runner = jackhmmer.Jackhmmer(
                binary_path=jackhmmer_binary_path,
                database_path=uniref90_database_path,
                get_tblout=True)

        if len(bfd_database_path) > 0:
            self.hhblits_bfd_runner = hhblits.HHBlits(
                binary_path=hhblits_binary_path,
                databases=[bfd_database_path])

        if len(uniref30_database_path) > 0:
            self.hhblits_uniref_runner = hhblits.HHBlits(
                binary_path=hhblits_binary_path,
                databases=[uniref30_database_path])

        if len(mgnify_database_path) > 0:
            self.jackhmmer_mgnify_runner = jackhmmer.Jackhmmer(
                binary_path=jackhmmer_binary_path,
                database_path=mgnify_database_path,
                get_tblout=True)

        if len(uniclust30_database_path) > 0:
            self.hhblits_uniclust_runner = hhblits.HHBlits(
                binary_path=hhblits_binary_path,
                databases=[uniclust30_database_path])

        if len(uniprot_database_path) > 0:
            self.jackhmmer_uniprot_runner = jackhmmer.Jackhmmer(
                binary_path=jackhmmer_binary_path,
                database_path=uniprot_database_path,
                get_tblout=True)

        # if os.path.exists(uniclust30_database_path):
        #     self.hhblits_uniclust_folddock_runner = hhblits.HHBlits(
        #         binary_path=hhblits_binary_path,
        #         databases=[uniclust30_database_path],
        #         all_seqs=True)

        if len(small_bfd_database_path) > 0:
            self.jackhmmer_small_bfd_runner = jackhmmer.Jackhmmer(
                binary_path=jackhmmer_binary_path,
                database_path=small_bfd_database_path,
                get_tblout=True)

        # if len(uniref30_database_path) > 0 and len(bfd_database_path) > 0:
        #     self.rosettafold_msa_runner = RosettaFold_Msa_runner(
        #         hhblits_binary_path=hhblits_binary_path,
        #         hhfilter_binary_path=hhfilter_binary_path,
        #         uniref30_database_path=uniref30_database_path,
        #         bfd_database_path=bfd_database_path)

        if len(uniclust30_database_path) > 0 and len(bfd_database_path) > 0:
            self.unclust30_bfd_msa_runner = hhblits.HHBlits(
                binary_path=hhblits_binary_path,
                databases=[bfd_database_path, uniclust30_database_path])

        if len(uniref30_database_path) > 0 and len(bfd_database_path) > 0:
            self.uniref30_bfd_msa_runner = hhblits.HHBlits(
                binary_path=hhblits_binary_path,
                databases=[bfd_database_path, uniref30_database_path])

        if len(colabfold_databases) > 0:
            self.colabfold_msa_runner = ColabFold_Msa_runner(colabfold_search_binary_path=colabfold_search_binary,
                                                             colabfold_split_msas_binary_path=colabfold_split_msas_binary,
                                                             mmseq_binary_path=mmseq_binary,
                                                             colabfold_databases=colabfold_databases)
        
        if len(JGIclust_database_path) > 0:
            self.deepmsa2_runner = DeepMSA2_runner(tool_path=deepmsa2_path,
                                                   bfd_database_path=bfd_database_path,
                                                   metaclust_database_path=metaclust_database_path,
                                                   mgnify_database_path=mgnify_database_path,
                                                   uniref90_database_path=uniref90_database_path,
                                                   uniref30_database_path=uniref30_database_path,
                                                   uniclust30_database_path=uniclust30_database_path,
                                                   JGIclust_database_path=JGIclust_database_path)

        if len(dhr_database_path) > 0:
            self.dhr_runner = DHR_runner(DHR_program_path=dhr_program_path, DHR_database_path=dhr_database_path)

    def process(self, input_fasta_path, msa_output_dir, multiprocess=True):
        """Runs alignment tools on the input sequence and creates features."""

        os.system(f"cp {input_fasta_path} {msa_output_dir}")
        
        targetname = pathlib.Path(input_fasta_path).stem

        msa_process_list = []

        if self.jackhmmer_uniref90_runner is not None:
            msa_process_list.append(
                [self.jackhmmer_uniref90_runner, input_fasta_path,
                 msa_output_dir, f'{targetname}_uniref90.sto', 'uniref90_sto'])

        if self.jackhmmer_mgnify_runner is not None:
            msa_process_list.append([self.jackhmmer_mgnify_runner, input_fasta_path,
                                     msa_output_dir, f'{targetname}_mgnify.sto', 'mgnify_sto'])

        if self.jackhmmer_small_bfd_runner is not None:
            msa_process_list.append(
                [self.jackhmmer_small_bfd_runner, input_fasta_path, msa_output_dir,
                 f'{targetname}_smallbfd.sto', 'smallbfd_sto'])

        if self.hhblits_bfd_runner is not None:
            msa_process_list.append([self.hhblits_bfd_runner, input_fasta_path,
                                     msa_output_dir, f'{targetname}_bfd.a3m', 'bfd_a3m'])

        if self.hhblits_uniref_runner is not None:
            msa_process_list.append([self.hhblits_uniref_runner, input_fasta_path,
                                     msa_output_dir, f'{targetname}_uniref30.a3m', 'uniref30_a3m'])

        if self.hhblits_uniclust_runner is not None:
            msa_process_list.append(
                [self.hhblits_uniclust_runner, input_fasta_path,
                 msa_output_dir, f'{targetname}_uniclust30.a3m', 'uniclust30_a3m'])

        if self.hhblits_uniclust_folddock_runner is not None:
            msa_process_list.append([self.hhblits_uniclust_folddock_runner, input_fasta_path,
                                    msa_output_dir, f'{targetname}_uniclust30_all.a3m', 'uniclust30_all_a3m'])

        if self.jackhmmer_uniprot_runner is not None:
            msa_process_list.append([self.jackhmmer_uniprot_runner, input_fasta_path,
                                    msa_output_dir, f'{targetname}_uniprot.sto', 'uniprot_sto'])

        if self.rosettafold_msa_runner is not None:
            msa_process_list.append(
                [self.rosettafold_msa_runner, input_fasta_path,
                msa_output_dir, f'{targetname}_rosettafold.a3m', 'rosettafold_sto'])

        if self.colabfold_msa_runner is not None:
            msa_process_list.append(
                [self.colabfold_msa_runner, input_fasta_path, msa_output_dir,
                 f'{targetname}_colabfold.a3m', 'colabfold_a3m'])

        if self.unclust30_bfd_msa_runner is not None:
            msa_process_list.append(
                [self.unclust30_bfd_msa_runner, input_fasta_path, msa_output_dir,
                 f'{targetname}_uniclust30_bfd.a3m', 'uniclust30_bfd_a3m'])

        if self.uniref30_bfd_msa_runner is not None:
            msa_process_list.append(
                [self.uniref30_bfd_msa_runner, input_fasta_path,
                 msa_output_dir, f'{targetname}_uniref30_bfd.a3m', 'uniref30_bfd_a3m'])

        if self.deepmsa2_runner is not None:
            msa_process_list.append(
                [self.deepmsa2_runner, input_fasta_path,
                 msa_output_dir, f'{targetname}_DeepMSA2.a3m', 'DeepMSA2_a3m'])

        if self.dhr_runner is not None:
            msa_process_list.append(
                [self.dhr_runner, input_fasta_path, msa_output_dir, 
                f'{targetname}_dhr.a3m', 'dhr_a3m'])

        if multiprocess:
            pool = Pool(processes=len(msa_process_list))
            results = pool.map(run_msa_tool, msa_process_list)
            pool.close()
            pool.join()
            #
            result_dict = {}
            for result in results:
                msa_key, msa_out_path = result
                if os.path.exists(msa_out_path):
                    result_dict[msa_key] = msa_out_path
        else:
            result_dict = {}
            for msa_process_params in msa_process_list:
                msa_key, msa_out_path = run_msa_tool(msa_process_params)
                if os.path.exists(msa_out_path):
                    result_dict[msa_key] = msa_out_path

        # Need to combine dhr a3m with other default a3ms
        dhr_af_a3m = os.path.join(msa_output_dir, f'{targetname}_dhr_af.a3m')
        if not os.path.exists(dhr_af_a3m):
            with open(result_dict['uniref90_sto']) as f:
                uniref90_msa = parsers.parse_stockholm(f.read())
                uniref90_msa = uniref90_msa.truncate(max_seqs=self.uniref_max_hits)

            with open(result_dict['mgnify_sto']) as f:
                mgnify_msa = parsers.parse_stockholm(f.read())
                mgnify_msa = mgnify_msa.truncate(max_seqs=self.mgnify_max_hits)

            with open(result_dict['uniref30_bfd_a3m']) as f:
                bfd_msa = parsers.parse_a3m(f.read())

            with open(result_dict['dhr_a3m']) as f:
                dhr_msa = parsers.parse_a3m(f.read())

            if len(dhr_msa.sequences) > 100000:
                cmd = f"{self.hhfilter_binary_path} -diff 50000 -i {result_dict['dhr_a3m']} -o {result_dict['dhr_a3m']}.filt -id 90"
                print(cmd)
                os.system(cmd)
                with open(result_dict['dhr_a3m'] + '.filt') as f:
                    dhr_msa = parsers.parse_a3m(f.read())

            combine_a3ms([dhr_msa, uniref90_msa, bfd_msa, mgnify_msa], dhr_af_a3m)

        result_dict['dhr_af_a3m'] = dhr_af_a3m

        return result_dict


class Monomer_alignment_generation_pipeline_img:
    """Runs the alignment tools and assembles the input features."""

    def __init__(self,
                 deepmsa_binary_path,
                 bfd_database_path,
                 img_database_path,
                 metaclust_database_path,
                 mgnify_database_path,
                 uniref90_database_path):
        """Initializes the data pipeline."""

        # alignment generation pipeline from alphafold
        self.img_msa_runner = IMG_Msa_runner(binary_path=deepmsa_binary_path,
                                             bfd_database_path=bfd_database_path,
                                             img_database_path=img_database_path,
                                             metaclust_database_path=metaclust_database_path,
                                             mgnify_database_path=mgnify_database_path,
                                             uniref90_database_path=uniref90_database_path)

    def process(self, input_fasta_path, msa_output_dir):
        """Runs alignment tools on the input sequence and creates features."""

        targetname = open(input_fasta_path).readlines()[0].rstrip('\n').lstrip('>')

        img_out_path = os.path.join(msa_output_dir, f'{targetname}.a3m')
        print(img_out_path)

        if not os.path.exists(img_out_path) or len(open(img_out_path).readlines()) == 0:
            img_out_path = self.img_msa_runner.query(input_fasta_path, msa_output_dir)

        return img_out_path
