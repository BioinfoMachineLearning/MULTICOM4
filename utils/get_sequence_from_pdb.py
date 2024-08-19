import os, sys, argparse, time


def get_sequence(inpdb):
    """Enclosing logic in a function to simplify code"""

    chain_sequences = {}
    res_codes = [
        # 20 canonical amino acids
        ('CYS', 'C'), ('ASP', 'D'), ('SER', 'S'), ('GLN', 'Q'),
        ('LYS', 'K'), ('ILE', 'I'), ('PRO', 'P'), ('THR', 'T'),
        ('PHE', 'F'), ('ASN', 'N'), ('GLY', 'G'), ('HIS', 'H'),
        ('LEU', 'L'), ('ARG', 'R'), ('TRP', 'W'), ('ALA', 'A'),
        ('VAL', 'V'), ('GLU', 'E'), ('TYR', 'Y'), ('MET', 'M'),
        # Non-canonical amino acids
        # ('MSE', 'M'), ('SOC', 'C'),
        # Canonical xNA
        ('  U', 'U'), ('  A', 'A'), ('  G', 'G'), ('  C', 'C'),
        ('  T', 'T'),
    ]

    three_to_one = dict(res_codes)
    # _records = set(['ATOM  ', 'HETATM'])
    _records = set(['ATOM  '])

    chain_max_res_num = {}
    for line in open(inpdb):
        line = line.strip()
        if line[0:6] in _records:
            resn = line[17:20]
            resi = int(line[22:26])
            icode = line[26]
            chain_id = line[21]
            if chain_id not in chain_max_res_num:
                chain_max_res_num[chain_id] = 0
            
            if resi > chain_max_res_num[chain_id]:
                chain_max_res_num[chain_id] = resi

    chain_sequences = {chain_id: ['-'] * chain_max_res_num[chain_id] for chain_id in chain_max_res_num}
    read = set()
    for line in open(inpdb):
        line = line.strip()
        if line[0:6] in _records:
            resn = line[17:20]
            resi = int(line[22:26])
            icode = line[26]
            chain_id = line[21]
            aa_resn = three_to_one.get(resn, 'X')
            chain_sequences[chain_id][resi-1] = aa_resn

    return chain_sequences


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--inpdb', type=str, required=True)

    args = parser.parse_args()

    chain_sequences = get_sequence(args.inpdb)
    for chain_id in sorted(chain_sequences.keys()):
        print(''.join(chain_sequences[chain_id]))
