# Modified from Alphafold2 codes

"""Library to run HHblits from Python."""

import glob
import os
import subprocess
from typing import Any, Mapping, Optional, Sequence
from absl import logging
from multicom4.tool import utils


class DHR_runner:
    """Python wrapper of the HHblits binary."""

    def __init__(self,
                 DHR_program_path,
                 DHR_database_path):

        self.DHR_program_path = DHR_program_path
        self.DHR_database_path = DHR_database_path

        print(f"Using database: {DHR_database_path}")

        if not os.path.exists(DHR_database_path):
            logging.error('Could not find DHR database %s', DHR_database_path)
            raise ValueError('Could not find DHR database %s', DHR_database_path)


    def query(self, input_fasta_path: str, output_a3m_path: str) -> Mapping[str, Any]:
        """Queries the database using DHR."""

        # targetname = open(input_fasta_path).readlines()[0].rstrip('\n').lstrip('>')

        # outpath = os.path.dirname(os.path.abspath(output_a3m_path))

        if not os.path.exists(output_a3m_path):
            cmd = [self.DHR_binary_path, self.DHR_program_path,
                '--input_path', input_fasta_path,
                '--database_path', self.DHR_database_path,
                '--output_a3m_path', output_a3m_path
            ]

            logging.info('DHR subprocess "%s"', ' '.join(cmd))

            process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

            with utils.timing('DHR query'):
                stdout, stderr = process.communicate()
                retcode = process.wait()

            if retcode:
                # Logs have a 15k character limit, so log HHblits error line by line.
                logging.error('DHR failed. DHR stderr begin:')
                for error_line in stderr.decode('utf-8').splitlines():
                    if error_line.strip():
                        logging.error(error_line.strip())
                logging.error('DHR stderr end')
                raise RuntimeError('DHR failed\nstdout:\n%s\n\nstderr:\n%s\n' % (
                    stdout.decode('utf-8'), stderr[:500_000].decode('utf-8')))

        return dict(a3m=output_a3m_path)
