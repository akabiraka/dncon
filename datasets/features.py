import sys
sys.path.append("../dncon")
import numpy as np
import pandas as pd

from datasets.atchley_factors import AtchleyFactors
from datasets.pssm import PSSM
from datasets.scratch import SCRATCH
from datasets.global_features import GlobalFeatures

class Features(object):
    def __init__(self, window_pos=5):
        super(Features, self).__init__()
        self.window_pos = window_pos
        self.atchley_factors = AtchleyFactors()
        self.pssm = PSSM()
        self.scratch = SCRATCH()
        self.global_features = GlobalFeatures()

    def compute(self, pdb_code, chain_id):
        pssm_df = self.pssm.get_pssm_observed_information(pdb_code, chain_id)
        ss3_df = self.scratch.get_secondary_structure_df(pdb_code, chain_id)
        acc2_df = self.scratch.get_solvency_accessibility_df(pdb_code, chain_id)
        exposed_prcnt_df = self.scratch.get_exposed_percent_df(pdb_code, chain_id)
        encoded_prot_len = self.global_features.encode_protein_length(pssm_df.shape[0])
        encoded_window_pos = self.global_features.encode_window_position(self.window_pos, pssm_df.shape[0])
        result = pd.concat([pssm_df, ss3_df, acc2_df, exposed_prcnt_df, encoded_prot_len, encoded_window_pos], axis=1, ignore_index=True)
        result = result.set_index(0)

        factors_df = pd.DataFrame()
        for aa in pssm_df.loc[:, 1]:
            factors_df = factors_df.append(self.atchley_factors.get(aa))
        
        result = pd.concat([result, factors_df], axis=1, ignore_index=True)       
        print("PSSM: {}\nSS3: {}\nACC2: {}\nExposed %: {}\nProtein len: {}\nWindow-pos: {}\nAtchley factors: {}\nCombined features: {}\n "\
            .format(pssm_df.shape, ss3_df.shape, acc2_df.shape, exposed_prcnt_df.shape, \
                encoded_prot_len.shape, encoded_window_pos.shape, factors_df.shape, result.shape))
        
        return result



features = Features()
features.compute("5sy8", "O")

