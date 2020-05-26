import sys
sys.path.append("../dncon")
import numpy as np
import pandas as pd
import configs.general_config as CONFIGS

class PSSM(object):
    def __init__(self):
        super(PSSM, self).__init__()

    def compute_things(self, pdb_code, chain_id):
        path = CONFIGS.PSIBLAST_DIR + pdb_code + chain_id + CONFIGS.DOT_PSSM
        df = pd.read_csv(path, delim_whitespace=True, header=None, skiprows=[1, 2]) # skipping 1st 2 rows
        df = df.head(-5) # removing last 5 rows
        return df
        
    def get_pssm(self, pdb_code, chain_id):
        """
        from psiblast output
        position-specific scoring matrix computed
        """
        df = self.compute_things(pdb_code, chain_id)
        pssm = df.loc[:, 2:21]
        return pssm

    def get_observed_percentage(self, pdb_code, chain_id):
        """
        from psiblast output
        weighted observed percentages rounded down
        """
        df = self.compute_things(pdb_code, chain_id)
        observed_percentage = df.loc[:, 22:41]
        return observed_percentage
    
    def get_information(self, pdb_code, chain_id):
        """
        from psiblast output
        information per position
        """
        df = self.compute_things(pdb_code, chain_id)
        information = df.loc[:, 42]
        return information

    def get_gapless_match(self, pdb_code, chain_id):
        """
        from psiblast output
        relative weight of gapless real matches to pseudocounts
        """
        df = self.compute_things(pdb_code, chain_id)
        gapless_match = df.loc[:, 43]
        return gapless_match
    
    def get_full(self, pdb_code, chain_id):
        """
        returns full dataframe
        """
        df = self.compute_things(pdb_code, chain_id)
        full = df.loc[:, 2:]
        return full

    def get_pssm_observed_percentage_information(self, pdb_code, chain_id):
        """
        returns pssm, weighted observed percentages and
        information per position
        """
        df = self.compute_things(pdb_code, chain_id)
        full_minus_gapless_match = df.loc[:, 2:42]
        return full_minus_gapless_match 



pssm = PSSM()
print(pssm.get_pssm_observed_percentage_information("5sy8", "O"))
print(pssm.get_pssm_observed_percentage_information("6y2d", "A"))


    