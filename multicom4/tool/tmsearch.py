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

"""Library to run HHsearch from Python."""

import glob
import os
import subprocess
from typing import Any, Dict, Mapping, Optional, Sequence, Tuple

from absl import logging
from multicom4.tool import utils
import pandas as pd
import pathlib
import datetime

# Internal import (7716).

class TMSearch:
    """Python wrapper of the TMsearch binary."""

    def __init__(self,
                 *,
                 binary_path: str,
                 program_path: str,
                 tmsearch_database: str):

        self.binary_path = binary_path
        self.program_path = program_path
        self.tmsearch_database = tmsearch_database

        if not os.path.exists(self.tmsearch_database):
            logging.error('Could not find PDB database %s for TMSearch', tmsearch_database)
            raise ValueError(f'Could not find PDB database {tmsearch_database} for TMSearch')
        
        self.PDBrepresent = os.path.join(program_path, 'representive_list')
        self.PDBallcluster = os.path.join(program_path, 'PDB_all_cluster')
        self.PDBcluster = os.path.join(program_path, 'PDB_noredundant_cluster')

    def query(self, pdb: str, outdir: str) -> str:

        if not os.path.exists(pdb):
            logging.error('Could not find input pdb: %s', pdb)
            raise ValueError(f'Could not find input pdb: {pdb}')

        tempfile = os.path.join(outdir, 'tmsearch.result')

        if not os.path.exists(tempfile):

            cmd = [self.binary_path, 
                pdb,
                '-dirc', self.tmsearch_database,
                self.PDBcluster, 
                '-outfmt', '1',
                ]

            logging.info('Launching subprocess "%s"', ' '.join(cmd))

            with open(tempfile, 'w') as fw:
                process = subprocess.Popen(cmd, stdout=fw, stderr=subprocess.PIPE)
                with utils.timing('TMSearch query'):
                    stdout, stderr = process.communicate()
                    retcode = process.wait()

            if retcode:
                # Stderr is truncated to prevent proto size errors in Beam.
                raise RuntimeError(
                    'HHSearch failed:\nstdout:\n%s\n\nstderr:\n%s\n' % (
                        stdout.decode('utf-8'), stderr[:100_000].decode('utf-8')))

        print("1111111111111111111111")

        data_dict = {'query': [], 'target': [], 'qaln': [], 'taln': [], 'qstart': [], 'qend': [], 'tstart': [], 'tend': [], 'evalue': [], 'alnlen': []}

        template = ''

        contents = open(tempfile).readlines()
        index = 0
        while index + 6 < len(contents):
            pdbline = contents[index]
            qaln = contents[index+1]
            templateline = contents[index+2]
            taln = contents[index+3]
            aln_detail = contents[index+4]
            _, _ = contents[index+5], contents[index+6]
            index += 7

            data_dict['query'] += [os.path.basename(pdb)]

            data_dict['target'] += [os.path.basename(templateline.split()[0])]
            data_dict['qaln'] += [qaln.rstrip('\n')]
            data_dict['taln'] += [taln.rstrip('\n')]

            data_dict['qstart'] += [0]
            sequence = [char for char in qaln.rstrip('\n') if char != '-']
            data_dict['qend'] += [len(sequence)]

            sequence = [char for char in taln.rstrip('\n') if char != '-']
            data_dict['tstart'] += [0]
            data_dict['tend'] += [len(sequence)]

            data_dict['evalue'] += [float(pdbline.rstrip('\n').split('TM-score=')[1])]
            data_dict['alnlen'] += [int(aln_detail.split()[1].split('Lali=')[1])]

        # outfile = os.path.join(outdir, 'tmsearch.csv')
        # pd.DataFrame(data_dict).to_csv(outfile)
        
        return pd.DataFrame(data_dict)