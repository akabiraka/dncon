import sys
sys.path.append("../dncon")
import numpy as np
import pandas as pd
import configs.general_config as CONFIGS

class PSSM(object):
    def __init__(self):
        super().__init__()

    def compute_things(self, pdb_code, chain_id):
        path = CONFIGS.PSIBLAST_DIR + pdb_code + chain_id + CONFIGS.DOT_PSSM
        df = pd.read_csv(path, delim_whitespace=True, header=None, skiprows=[1, 2]) # skipping 1st 2 rows
        df = df.head(-5) # removing last 5 rows
        self.pssm = df.loc[:, 2:21]
        self.observed_percentage = df.loc[:, 22:41]
        self.information = df.loc[:, 42]
        self.gapless_match = df.loc[:, 43]
        self.full = df.loc[:, 2:]
        self.full_minus_gapless_match = df.loc[:, 2:42]
        
    def get_pssm(self, pdb_code, chain_id):
        """
        from psiblast output
        position-specific scoring matrix computed
        """
        self.compute_things(pdb_code, chain_id)
        return self.pssm

    def get_observed_percentage(self, pdb_code, chain_id):
        """
        from psiblast output
        weighted observed percentages rounded down
        """
        self.compute_things(pdb_code, chain_id)
        return self.observed_percentage
    
    def get_information(self, pdb_code, chain_id):
        """
        from psiblast output
        information per position
        """
        self.compute_things(pdb_code, chain_id)
        return self.information

    def get_gapless_match(self, pdb_code, chain_id):
        """
        from psiblast output
        relative weight of gapless real matches to pseudocounts
        """
        self.compute_things(pdb_code, chain_id)
        return self.gapless_match
    
    def get_full(self, pdb_code, chain_id):
        """
        returns full dataframe
        """
        self.compute_things(pdb_code, chain_id)
        return self.full

    def get_full_minus_last_col(self, pdb_code, chain_id):
        """
        returns pssm, weighted observed percentages and
        information per position
        """
        self.compute_things(pdb_code, chain_id)
        return self.full_minus_gapless_match 



pssm = PSSM()
print(pssm.get_pssm("5sy8", "O"))
print(pssm.get_pssm("6y2d", "A"))


    