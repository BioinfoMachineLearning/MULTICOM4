import os, sys, argparse, time
from multiprocessing import Pool
from tqdm import tqdm
import pandas as pd
import numpy as np
from multicom4.monomer_templates_concatenation import parsers
from multicom4.monomer_alignment_generation.alignment import *
from multicom4.common.pipeline import run_monomer_msa_pipeline

mgnify_max_hits = 501
uniref_max_hits = 10000


def retrieve_sequence_ids_domain(ids, regex=None):
    if regex is None:
        regex = ID_EXTRACTION_REGEX

    sequence_ids = []
    headers = []
    id_to_full_header = defaultdict(list)

    for current_id in ids:
        match = False
        for pattern in regex:
            m = re.search(pattern, current_id)
            # require a non-None match and at least one extracted pattern
            if m and len(m.groups()) > 0:
                # this extracts the parenthesized match
                sequence_ids.append(m.group(1))
                id_to_full_header[m.group(1)].append(current_id)
                headers += [current_id]
                match = True
                break
        
        if not match:
            header = current_id.split()[0].split('/')[0]
            sequence_ids.append(header)
            id_to_full_header[header].append(current_id)
            headers += [current_id]

    return sequence_ids, id_to_full_header, headers

class DomainAlignment:
    def __init__(self, ids, seqs, annotation=None):

        if len(ids) == 0:
            raise ValueError("Alignment file is empty!")

        self.main_id = ids[0]
        self.main_seq = seqs[0]
        self.L = len(seqs[0])

        self.ids, self.headers_dict, self.headers = retrieve_sequence_ids_domain(ids[1:])
        if len(self.ids) == 0:
            self.ids = ids[1:]
            self.headers = defaultdict(list)
            for id in self.ids:
                self.headers_dict[id].append(id)
                self.headers += [id]
                
        self.seqs = seqs[1:]
        self.N = len(self.seqs)

        if len(self.ids) != self.N:
            print(self.ids)
            print(ids)
            raise ValueError(
                f"Number of sequence IDs and length of alignment do not match: {self.N} and {len(self.ids)}"
            )
        self.id_to_index = {id_: i for i, id_ in enumerate(self.ids)}

    @classmethod
    def from_dict(cls, seq_dict, annotation=None):
        return cls(list(seq_dict.keys()), list(seq_dict.values()), annotation)

    @classmethod
    def from_a3m(cls, ina3m):
        with open(ina3m) as fileobj:
            seqs = read_a3m(fileobj)
        return cls.from_dict(seqs)

    def __getitem__(self, index):
        """
        .. todo::
            eventually this should allow fancy indexing and offer the functionality of select()
        """
        if index in self.id_to_index:
            return self.seqs[self.id_to_index[index]]
        elif index in range(self.N):
            return self.seqs[index]
        else:
            raise KeyError(
                "Not a valid index for sequence alignment: {}".format(index)
            )

    def __len__(self):
        return self.N

def make_msa_features(msas, msa_save_path):
    """Constructs a feature dict of MSA features."""
    if not msas:
        raise ValueError('At least one MSA must be provided.')

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

    with open(msa_save_path, 'w') as fw:
        for (desc, seq) in zip(seen_desc, seen_sequences):
            fw.write(f'>{desc}\n{seq}\n')


def _make_msa_df(domain_alignment):
    msa_df = pd.DataFrame({
        'msa_ids': list(domain_alignment.ids),
        'msa_row': np.arange(len(list(domain_alignment.ids))),
    })
    return msa_df


def _create_id_dict(msa_df):
    """Creates mapping from species to msa dataframe of that species."""
    species_lookup = {}
    for species, species_df in msa_df.groupby('msa_ids'):
        species_lookup[species] = species_df
    return species_lookup


def _match_rows_by_seqid(this_id_msa_dfs):
    all_paired_msa_rows = []
    num_seqs = [len(id_df) for id_df in this_id_msa_dfs if id_df is not None]
    take_num_seqs = np.min(num_seqs)
    for id_df in this_id_msa_dfs:
        if id_df is not None:
            msa_rows = id_df.msa_row.iloc[:take_num_seqs].values
        else:
            msa_rows = [-1] * take_num_seqs  # take the last 'padding' row
        all_paired_msa_rows.append(msa_rows)
    all_paired_msa_rows = list(np.array(all_paired_msa_rows).transpose())
    return all_paired_msa_rows

def reorder_paired_rows(all_paired_msa_rows_dict) :
    all_paired_msa_rows = []
    for num_pairings in sorted(all_paired_msa_rows_dict, reverse=True):
        paired_rows = all_paired_msa_rows_dict[num_pairings]
        paired_rows_product = abs(np.array([np.prod(rows) for rows in paired_rows]))
        paired_rows_sort_index = np.argsort(paired_rows_product)
        all_paired_msa_rows.extend(paired_rows[paired_rows_sort_index])
    return np.array(all_paired_msa_rows)

def write_paired_a3ms(sequence, domain_alignments, domain_ranges, paired_rows, outfile):
    def _prepare_header(ids):
        return "_____".join(ids)

    sequences_to_write = {}
    filter_pair_ids = {}
    for i in range(len(domain_alignments)):
        filter_pair_ids[f"id_{i + 1}"] = [domain_alignments[i].main_id]
        filter_pair_ids[f"index_{i + 1}"] = [0]

    pair_id_count = 0
    for pair_index in range(1, len(list(paired_rows[:, 0]))):

        row_indices = list(paired_rows[pair_index, :])
        seqs = []
        headers = []

        for j in range(len(domain_alignments)):
            index = row_indices[j]
            header = ""
            if index == -1:
                header = f'placeholder{pair_id_count}'
                seq = '-' * len(domain_alignments[j].main_seq)
            else:
                seq = domain_alignments[j].seqs[index]
                header = domain_alignments[j].ids[index]

            headers += [header]
            seqs += [seq]
            filter_pair_ids[f'id_{j + 1}'] += [header]
            filter_pair_ids[f'index_{j + 1}'] += [pair_id_count + 1]

        concatenated_header = _prepare_header(headers)

        # save the information
        sequences_to_write[concatenated_header] = seqs

        pair_id_count += 1

    # add paired sequences
    sequences_full = {}
    for header_full in sequences_to_write:
        paired_sequence = ["-"] * len(sequence)
        for domain_sequence, (domain_starts, domain_ends) in zip(sequences_to_write[header_full], domain_ranges):
            domain_idx = 0
            domain_sequence_list = list(domain_sequence)
            for domain_start, domain_end in zip(domain_starts, domain_ends):
                print(domain_start)
                print(domain_end)
                domain_sequence_length = domain_end - domain_start + 1
                paired_sequence[domain_start:domain_end+1] = domain_sequence_list[domain_idx:domain_idx+domain_sequence_length]
                domain_idx += domain_sequence_length

        seq_full = "".join(paired_sequence)
        sequences_full[header_full] = seq_full

    # add unpaired sequences
    for domain_alignment, (domain_starts, domain_ends) in zip(domain_alignments, domain_ranges):
        #print(domain_alignment.ids)
        for seqid in domain_alignment.ids:
            domain_sequence = domain_alignment[seqid]
            unpaired_sequence = ["-"] * len(sequence)
            domain_idx = 0
            domain_sequence_list = list(domain_sequence)
            for domain_start, domain_end in zip(domain_starts, domain_ends):
                domain_sequence_length = domain_end - domain_start + 1
                unpaired_sequence[domain_start:domain_end+1] = domain_sequence_list[domain_idx:domain_idx+domain_sequence_length]
                domain_idx += domain_sequence_length
            seq_full = "".join(unpaired_sequence)
            sequences_full[f"{domain_alignment.main_id}_{seqid}"] = seq_full
            # print(sequences_full.keys())

    with open(outfile, "w") as of:
        write_a3m(sequences_full, of)

def run_hhsearch_dom(disorder_pred, N1_outdir, fasta_path, params):
    # 1. hhsearch
    print(f"Predicting domain boundries using hhsearch")
    dom_hhsearch_out = os.path.join(N1_outdir, 'Dom_by_hhsearch')
    domain_info_file = os.path.join(dom_hhsearch_out, 'domain_info')
    if os.path.exists(domain_info_file):
        print(f"Found {domain_info_file}!!")
    else:
        cmd = f"perl {params['dom_hhsearch_script']} {fasta_path} {dom_hhsearch_out} {params['hhsearch15db']} 40"
        print(cmd)
        os.system(cmd)
        os.system(f"cp {dom_hhsearch_out}/domain_info_thre40 {domain_info_file}")

    domain_disorder_info, domain_range_info = correct_domain(domain_info_file, disorder_pred)
    domain_def = os.path.join(dom_hhsearch_out, 'domain_info_final')
    with open(domain_def, 'w') as fw:
        print("Corrected domain regions:")
        for i in range(len(domain_range_info)):
            print(f"domain {i}: {domain_range_info[i]} {domain_disorder_info[i]}")
            fw.write(f"domain {i}: {domain_range_info[i]} {domain_disorder_info[i]}\n")

    return domain_def

def run_dom_parser(disorder_pred, N1_outdir, fasta_path, ranked_0_pdb, params):
    dom_parse_out = os.path.join(N1_outdir, 'Dom_by_domain_parser')
    domain_info_file = os.path.join(dom_parse_out, 'domain_info')
    print(f"Found input pdb, using domain parser to predict domain boundries")
    if len(sequence) < 3000:     
        os.makedirs(dom_parse_out, exist_ok=True)
        if os.path.exists(domain_info_file):
            print(f"Found {domain_info_file}!!")
        else:
            cmd = f"perl {params['dom_parser_script']} {fasta_path} {ranked_0_pdb} {dom_parse_out}"
            print(cmd)
            os.system(cmd)
    else:
        print(f"Domain parser cannot deal with sequence longer than 3000, skip!")
        return None

    domain_disorder_info, domain_range_info = correct_domain(domain_info_file, disorder_pred)
    domain_def = os.path.join(dom_parse_out, 'domain_info_final')
    with open(domain_def, 'w') as fw:
        print("Corrected domain regions:")
        for i in range(len(domain_range_info)):
            print(f"domain {i}: {domain_range_info[i]} {domain_disorder_info[i]}")
            fw.write(f"domain {i}: {domain_range_info[i]} {domain_disorder_info[i]}\n")

    return domain_def


def run_unidoc(disorder_pred, N1_outdir, fasta_path, ranked_0_pdb, params):
    unidoc_parse_out = os.path.join(N1_outdir, 'Dom_by_UniDoc')
    domain_info_file = os.path.join(unidoc_parse_out, 'domain_info')
    print(f"Found input pdb, using UniDoc to predict domain boundries")
    os.makedirs(unidoc_parse_out, exist_ok=True)
    if os.path.exists(domain_info_file):
        print(f"Found {domain_info_file}!!")
    else:
        os.chdir(params['unidoc_program_path'])
        cmd = f"python {params['unidoc_program_script']} -i {ranked_0_pdb} -c A > {unidoc_parse_out}/unidoc.out"
        print(cmd)
        os.system(cmd)

        domain_range_out = open(os.path.join(unidoc_parse_out, 'unidoc.out')).readlines()[0]

        domain_ranges = domain_range_out.rstrip('\n').split('/')

        sorted_domain_ranges = sorted(domain_ranges, key=lambda x: int(x.split('~')[0]))
        print(sorted_domain_ranges)

        with open(domain_info_file, 'w') as fw:
            for i, domain_str in enumerate(sorted_domain_ranges):
                domain_str = domain_str.replace('~', '-')
                fw.write(f"domain {i}: {domain_str}\n")

    domain_disorder_info, domain_range_info = correct_domain(domain_info_file, disorder_pred)
    domain_def = os.path.join(unidoc_parse_out, 'domain_info_final')
    with open(domain_def, 'w') as fw:
        print("Corrected domain regions:")
        for i in range(len(domain_range_info)):
            print(f"domain {i}: {domain_range_info[i]} {domain_disorder_info[i]}")
            fw.write(f"domain {i}: {domain_range_info[i]} {domain_disorder_info[i]}\n")

    return domain_def


def generate_domain_fastas(fasta_path, domain_info, output_dir):

    os.makedirs(outputdir, exist_ok=True)

    sequence = open(fasta_path).readlines()[1].rstrip('\n')
    for line in open(domain_info):
        line = line.rstrip('\n')
        # domain 1: 140-317 Normal
        # domain 1: 140-317
        domain_name, domain_range = line.split(':')
        domain_name = domain_name.replace(' ', '')
        domain_range = domain_range.split()[0]

        domain_sequence = ""
        starts, ends = [], []

        for domain_range_str in domain_range.split(','):
            start, end = domain_range_str.split('-')
            starts += [int(start)-1]
            ends += [int(end)]

        for start, end in zip(starts, ends):
            domain_sequence += sequence[start:end]

        with open(os.path.join(outputdir, domain_name + '.fasta'), 'w') as fw:
            fw.write(f">{domain_name}\n{domain_sequence}")


def generate_domain_alignments(params, domain_fastas_dir, output_dir):

    makedir_if_not_exists(output_dir)

    for fasta_path in os.listdir(domain_fastas_dir):

        fasta_path = os.path.join(domain_fastas_dir, fasta_path)

        # targetname = pathlib.Path(FLAGS.fasta_path).stem
        targetname = None
        sequence = None
        for line in open(fasta_path):
            line = line.rstrip('\n').strip()
            if line.startswith('>'):
                targetname = line[1:].split()[0]
                # if targetname_in_fasta != targetname:
                #     print("Warning: fasta file name doesn't match with fasta content!")
            else:
                sequence = line

        outdir = os.path.join(output_dir, targetname)

        makedir_if_not_exists(outdir)

        print("#################################################################################################")
        print(f"1. Start to generate alignments for monomers")

        uniref30 = params['uniref_db']
        uniclust30 = ''
        uniref90_fasta = params['uniref90_fasta']

        smallbfd = ""  # params['smallbfd_database']
        bfd = params['bfd_database']
        mgnify = params['mgnify_database']

        hhblits_binary = params['hhblits_program']
        hhfilter_binary = params['hhfilter_program']
        jackhmmer_binary = params['jackhmmer_program']

        colabfold_search_binary = '' #params['colabfold_search_program']
        colabfold_split_msas_binary = '' #params['colabfold_split_msas_program']
        colabfold_databases = '' #params['colabfold_databases']
        mmseq_binary = '' #params['mmseq_program']

        deepmsa2_path = '' #params['deepmsa2_path']
        JGIclust_database_path = '' #params['JGIclust_database']
        metaclust_database_path = '' #params['metaclust_database']
        
        dhr_program_path = '' #params['DHR_program_path']
        dhr_database_path = '' #params['DHR_database_path']

        result = run_monomer_msa_pipeline(fasta=fasta_path, outdir=outdir, params=params, only_monomer=True)

        if result is None:
            raise RuntimeError('The monomer alignment generation has failed!')

def combine_domain_msas(sequence, domain_info, domaindir):

    all_domain_id_dict = []
    all_domain_ranges = []
    domain_alignments = []
    common_ids = set()
    for line in open(domain_info):
        line = line.rstrip('\n')
        # domain 1: 140-317 Normal
        # domain 1: 140-317
        domain_name, domain_range = line.split(':')
        domain_name = domain_name.replace(' ', '')
        domain_range = domain_range.split()[0]

        domain_sequence = ""
        starts, ends = [], []

        for domain_range_str in domain_range.split(','):
            start, end = domain_range_str.split('-')

            # list index start from 0
            starts += [int(start) - 1]
            ends += [int(end) - 1]

        all_domain_ranges += [(starts, ends)]

        alndir = os.path.join(domaindir, domain_name)

        uniref_sto = os.path.join(alndir, domain_name + '_uniref90.sto')
        with open(uniref_sto) as f:
            uniref90_msa = parsers.parse_stockholm(f.read())
            uniref90_msa = uniref90_msa.truncate(max_seqs=uniref_max_hits)

        bfd_uniref30_a3m = os.path.join(alndir, domain_name + '_uniref30_bfd.a3m')
        with open(bfd_uniref30_a3m) as f:
            bfd_msa = parsers.parse_a3m(f.read())

        mgnify_a3m = os.path.join(alndir, domain_name + '_mgnify.sto')
        with open(mgnify_a3m) as f:
            mgnify_msa = parsers.parse_stockholm(f.read())
            mgnify_msa = mgnify_msa.truncate(max_seqs=mgnify_max_hits)

        domain_a3m_dir = os.path.join(domaindir, 'domain_a3ms')
        os.makedirs(domain_a3m_dir, exist_ok=True)

        alphafold_a3m = os.path.join(domain_a3m_dir, domain_name + '_alphafold.a3m')
        make_msa_features((uniref90_msa, bfd_msa, mgnify_msa), msa_save_path=alphafold_a3m)

        domain_alignment = DomainAlignment.from_a3m(alphafold_a3m)
        domain_alignments += [domain_alignment]

        msa_df = _make_msa_df(domain_alignment)
        id_dict = _create_id_dict(msa_df)
        all_domain_id_dict.append(id_dict)
        common_ids.update(set(id_dict))

    common_ids = sorted(common_ids)
    print(common_ids)

    num_examples = len(all_domain_id_dict)
    all_paired_msa_rows = [np.zeros(num_examples, int)]
    all_paired_msa_rows_dict = {k: [] for k in range(num_examples)}
    all_paired_msa_rows_dict[num_examples] = [np.zeros(num_examples, int)]

    for seqid in common_ids:
        if not seqid:
            continue
        this_id_msa_dfs = []
        id_dfs_present = 0
        for id_dict in all_domain_id_dict:
            if seqid in id_dict:
                this_id_msa_dfs.append(id_dict[seqid])
                id_dfs_present += 1
            else:
                this_id_msa_dfs.append(None)

        # Skip species that are present in only one chain.
        if id_dfs_present <= 1:
            continue

        if np.any(
                np.array([len(id_df) for id_df in
                          this_id_msa_dfs if
                          isinstance(id_df, pd.DataFrame)]) > 600):
            continue

        paired_msa_rows = _match_rows_by_seqid(this_id_msa_dfs)
        all_paired_msa_rows.extend(paired_msa_rows)
        all_paired_msa_rows_dict[id_dfs_present].extend(paired_msa_rows)

    paired_chains_to_paired_row_indices = {
        num_examples: np.array(paired_msa_rows) for
        num_examples, paired_msa_rows in all_paired_msa_rows_dict.items()
    }

    paired_rows = reorder_paired_rows(paired_chains_to_paired_row_indices)

    # print(paired_rows)

    outfile = os.path.join(outdir, 'paired.a3m')
    write_paired_a3ms(sequence=sequence, domain_alignments=domain_alignments, 
                      domain_ranges=all_domain_ranges, paired_rows=paired_rows, outfile=outfile)

    return outfile