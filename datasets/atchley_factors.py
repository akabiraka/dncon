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
        self.atchley_df = atchley_df.set_index(0)
        # print(atchley_df.loc['A'])

    def get(self, amino_acid):
        result = self.atchley_df.loc[amino_acid]
        return result


# af = AtchleyFactors()
# print(af.get('A'))
# to run 
# python datasets/atchley_factors.py