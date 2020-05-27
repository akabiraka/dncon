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

    def __get_one_hot(self, path, dict):
        """
        Private method
        Given the .acc or .ss file path and respective dictionary of characters,
        this method returns the 1-hot encoding of that data.
        """
        file_content = open(path, "r")
        data = list(file_content)[1].rstrip("\n")
        numeric = [dict[ch] for ch in data]
        one_hot_tensor = F.one_hot(torch.tensor(numeric), num_classes=len(dict))
        return one_hot_tensor

    def get_secondary_structure(self, pdb_code, chain_id):
        """
        returns predicted secondary-structure in tensor
        """
        path = CONFIGS.SCRATCH_OUT + pdb_code + chain_id + CONFIGS.DOT_SS
        return self.__get_one_hot(path, CONFIGS.SS3)
    
    def get_secondary_structure_df(self, pdb_code, chain_id):
        """
        returns predicted secondary-structure in dataframe
        """
        path = CONFIGS.SCRATCH_OUT + pdb_code + chain_id + CONFIGS.DOT_SS
        ss3_tensor = self.__get_one_hot(path, CONFIGS.SS3)
        ss3_df = pd.DataFrame(ss3_tensor.numpy())
        return ss3_df

    def get_solvency_accessibility(self, pdb_code, chain_id):
        """
        returns solvency-accessibility in tensor
        """
        path = CONFIGS.SCRATCH_OUT + pdb_code + chain_id + CONFIGS.DOT_ACC
        return self.__get_one_hot(path, CONFIGS.ACC2)

    def get_solvency_accessibility_df(self, pdb_code, chain_id):
        """
        returns solvency-accessibility in dataframe
        """
        path = CONFIGS.SCRATCH_OUT + pdb_code + chain_id + CONFIGS.DOT_ACC
        acc2_tensor = self.__get_one_hot(path, CONFIGS.ACC2)
        acc2_df = pd.DataFrame(acc2_tensor.numpy())
        return acc2_df

    def get_exposed_percent(self, pdb_code, chain_id):
        """ 
        Compute the percentage of predicted exposed alpha, beta and coli residues
        from secondary structure.
        returns in numpy format
        """
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
        """ 
        Compute the percentage of predicted exposed alpha, beta and coli residues
        from secondary structure.
        returns in dataframe format
        """
        exposed_percent = self.get_exposed_percent(pdb_code, chain_id)
        return pd.DataFrame(exposed_percent)

    def __get_data(self, path):
        """
        read a file and returns the second line by removing any newline
        """
        file_content = open(path, "r")
        data = list(file_content)[1].rstrip("\n")
        return data

# scratch = SCRATCH()
# print(scratch.get_secondary_structure("6y2d", "A"))
# print(scratch.get_exposed_percent("4eiu", "A"))