import sys
sys.path.append("../dncon")
import numpy as np
import pandas as pd
import torch
import torch.nn.functional as F

import configs.general_config as CONFIGS

class SCRATCH(object):
    def __init__(self):
        super(SCRATCH, self).__init__()

    def get_one_hot(self, path, dict):
        file_content = open(path, "r")
        data = list(file_content)[1].rstrip("\n")
        numeric = [dict[ch] for ch in data]
        one_hot_tensor = F.one_hot(torch.tensor(numeric), num_classes=len(dict))
        return one_hot_tensor

    def get_secondary_structure(self, pdb_code, chain_id):
        path = CONFIGS.SCRATCH_OUT + pdb_code + chain_id + CONFIGS.DOT_SS
        return self.get_one_hot(path, CONFIGS.SS3)
    
    def get_secondary_structure_df(self, pdb_code, chain_id):
        path = CONFIGS.SCRATCH_OUT + pdb_code + chain_id + CONFIGS.DOT_SS
        ss3_tensor = self.get_one_hot(path, CONFIGS.SS3)
        ss3_df = pd.DataFrame(ss3_tensor.numpy())
        return ss3_df

    def get_solvency_accessibility(self, pdb_code, chain_id):
        path = CONFIGS.SCRATCH_OUT + pdb_code + chain_id + CONFIGS.DOT_ACC
        return self.get_one_hot(path, CONFIGS.ACC2)

    def get_solvency_accessibility_df(self, pdb_code, chain_id):
        path = CONFIGS.SCRATCH_OUT + pdb_code + chain_id + CONFIGS.DOT_ACC
        acc2_tensor = self.get_one_hot(path, CONFIGS.ACC2)
        acc2_df = pd.DataFrame(acc2_tensor.numpy())
        return acc2_df

    def get_exposed_percent(self, pdb_code, chain_id):
        acc2_path = CONFIGS.SCRATCH_OUT + pdb_code + chain_id + CONFIGS.DOT_ACC
        ss3_path = CONFIGS.SCRATCH_OUT + pdb_code + chain_id + CONFIGS.DOT_SS
        acc2_data = self.__get_data(acc2_path)
        ss3_data = self.__get_data(ss3_path)
        exposed = ""
        for i, ch in enumerate(acc2_data):
            if ch == 'e':
               exposed = exposed + ss3_data[i]
        n_total = len(exposed)
        perct_C = exposed.count('C') / n_total
        perct_E = exposed.count('E') / n_total
        perct_H = exposed.count('H') / n_total
        exposed_percent = np.zeros((len(acc2_data), 3))
        exposed_percent[:] = np.array([perct_C, perct_E, perct_H])
        # print(broadcasted)
        return exposed_percent
        
    def get_exposed_percent_df(self, pdb_code, chain_id):
        exposed_percent = self.get_exposed_percent(pdb_code, chain_id)
        return pd.DataFrame(exposed_percent)

    def __get_data(self, path):
        file_content = open(path, "r")
        data = list(file_content)[1].rstrip("\n")
        return data

# scratch = SCRATCH()
# print(scratch.get_secondary_structure("6y2d", "A"))
# print(scratch.get_exposed_percent("4eiu", "A"))