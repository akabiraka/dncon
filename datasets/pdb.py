import sys
sys.path.append('../dncon')
import numpy as np

from Bio.PDB import *
from Bio import SeqIO

import configs.general_config as CONFIGS

class PDB(object):
    def __init__(self):
        print("alhumdulillah")
        self.pdb_ids_file = CONFIGS.ALL_PDB_IDS

        self.pdbl = PDBList()
        self.parser = MMCIFParser(QUIET=True)

    def read_pdb_ids(self, line):
        """
            Assuming a line has 6 columns in a row.
            returns:
                pdb_code, chain_id: ('5SY8', 'O') in upper case
        """
        line = line.split()
        pdb_code = line[0][:-1]
        chain_id = line[0][-1]
        return pdb_code, chain_id

    def download(self, pdb_code):
        # downloading
        self.pdbl.retrieve_pdb_file(pdb_code, pdir=CONFIGS.PDB_DIR, file_format=CONFIGS.CIF)

pdb = PDB()
file_content = open(CONFIGS.ALL_PDB_IDS, "r")
for i, line in enumerate(file_content):
    pdb_code, chain_id = pdb.read_pdb_ids(line)
    print(pdb_code, chain_id)
    pdb.download(pdb_code)

