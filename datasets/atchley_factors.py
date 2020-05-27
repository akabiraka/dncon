import sys
sys.path.append('../dncon')
import numpy as np
import pandas as pd
import configs.general_config as CONFIGS

class AtchleyFactors(object):
    def __init__(self):
        super(AtchleyFactors, self).__init__()
        atchley_file = CONFIGS.ATCHLEY_FACTORS
        atchley_df = pd.read_csv(atchley_file, delim_whitespace=True, header=None)
        atchley_df = atchley_df.drop(0)
        self.atchley_df = atchley_df.set_index(0) # setting index to the amino-acid characters
        # print(atchley_df.loc['A'])

    def get(self, amino_acid):
        """
        returns atchley factors for given amino-acid character
        """
        result = self.atchley_df.loc[amino_acid]
        return result

    def get_all(self, amino_acids):
        factors_df = pd.DataFrame()
        for aa in amino_acids:
            factors_df = factors_df.append(self.get(aa))
        factors_df.reset_index(drop=True, inplace=True)
        return factors_df


# af = AtchleyFactors()
# print(af.get('A'))
# to run 
# python datasets/atchley_factors.py