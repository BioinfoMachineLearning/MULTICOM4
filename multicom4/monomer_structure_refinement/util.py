import copy
import os
import sys
import time, json
from multicom4.common.util import makedir_if_not_exists, check_dirs
import dataclasses


def combine_a3ms(infiles, outfile):
    descriptions = []
    seqs = []
    for infile in infiles:
        for line in open(infile):
            line = line.rstrip('\n')
            if line.startswith('>'):
                descriptions += [line]
            else:
                seqs += [line]

    with open(outfile, 'w') as fw:
        for (desc, seq) in zip(descriptions, seqs):
            fw.write(f"{desc}\n{seq}\n")


def convert_taln_seq_to_a3m(query_non_gaps, aln):
    for is_query_res_non_gap, sequence_res in zip(query_non_gaps, aln):
        if is_query_res_non_gap:
            yield sequence_res


def cal_tmscore(tmscore_program, inpdb, nativepdb, tmpdir):
    cwd = os.getcwd()

    makedir_if_not_exists(tmpdir)

    os.chdir(tmpdir)

    os.system(f"cp {inpdb} inpdb.pdb")
    os.system(f"cp {nativepdb} native.pdb")

    inpdb = "inpdb.pdb"
    nativepdb = "native.pdb"

    cmd = tmscore_program + ' ' + inpdb + ' ' + nativepdb + " | grep TM-score | awk '{print $3}' "
    print(cmd)
    tmscore_contents = os.popen(cmd).read().split('\n')
    tmscore = float(tmscore_contents[2].rstrip('\n'))
    cmd = tmscore_program + ' ' + inpdb + ' ' + nativepdb + " | grep GDT-score | awk '{print $3}' "
    tmscore_contents = os.popen(cmd).read().split('\n')
    gdtscore = float(tmscore_contents[0].rstrip('\n'))

    os.chdir(cwd)

    os.system("rm -rf " + tmpdir)

    return tmscore, gdtscore


def build_alignment_indices(sequence, start_index):
    indices_list = []
    counter = start_index
    for symbol in sequence:
        if symbol == '-':
            indices_list.append(-1)
        else:
            indices_list.append(counter)
            counter += 1
    return indices_list
