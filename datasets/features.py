import sys
sys.path.append("../dncon")
import numpy as np
import pandas as pd

from atchley_factors import AtchleyFactors
from pssm import PSSM
from scratch import SCRATCH

class Features(object):
    def __init__(self):
        super(Features, self).__init__()
        self.atchley = AtchleyFactors()
        self.pssm = PSSM()
        self.scratch = SCRATCH()

    def compute(self, pdb_code, chain_id):
        pssm_df = self.pssm.get_pssm_observed_information(pdb_code, chain_id)
        ss3_df = self.scratch.get_secondary_structure_df(pdb_code, chain_id)
        acc2_df = self.scratch.get_solvency_accessibility_df(pdb_code, chain_id)

        result = pd.concat([pssm_df, ss3_df, acc2_df], axis=1)
        print(pssm_df.shape, ss3_df.shape, acc2_df.shape, result.shape)
        # print(result)

features = Features()
features.compute("5sy8", "O")

