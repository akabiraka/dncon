import sys
sys.path.append('../dncon')
import numpy as np
import pandas as pd
import configs.general_config as CONFIGS

class PairwisePotentials(object):
    def __init__(self):
        super(PairwisePotentials, self).__init__()
        self.brauns_df = pd.read_csv(CONFIGS.BRAUNS).set_index('aa')
        self.jernigans_df = pd.read_csv(CONFIGS.JERNIGANS).set_index('aa')
        self.levitts_df = pd.read_csv(CONFIGS.LEVITTS).set_index('aa')

    def get_brauns_potentials(self, aa1, aa2):
        # print(self.brauns_df[aa1][aa2], self.brauns_df[aa2][aa1])
        return self.brauns_df[aa1][aa2]

    def get_jernigans_potentials(self, aa1, aa2):
        # print(self.jernigans_df[aa1][aa2], self.jernigans_df[aa2][aa1])
        return self.jernigans_df[aa1][aa2]

    def get_levitts_potentials(self, aa1, aa2):
        # print(self.levitts_df[aa1][aa2], self.levitts_df[aa2][aa1])
        return self.levitts_df[aa1][aa2]

# pp = PairwisePotentials()
# print(pp.get_brauns_potentials('G', 'A'))
# print(pp.get_jernigans_potentials('G', 'A'))
# print(pp.get_levitts_potentials('G', 'A'))