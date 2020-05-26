import sys
sys.path.append("../dncon")
import numpy as np
import torch
import torch.nn.functional as F

import configs.general_config as CONFIGS

class Scratch(object):
    def __init__(self):
        super(Scratch, self).__init__()

    def get_one_hot(self, path, dict):
        file_content = open(path, "r")
        data = list(file_content)[1].rstrip("\n")
        numeric = [dict[ch] for ch in data]
        one_hot_tensor = F.one_hot(torch.tensor(numeric),
                                       num_classes=len(dict))
        return one_hot_tensor

    def get_secondary_structure(self, pdb_code, chain_id):
        path = CONFIGS.SCRATCH_OUT + pdb_code + chain_id + CONFIGS.DOT_SS
        return self.get_one_hot(path, CONFIGS.SS3)

    def get_solvency_accessibility(self, pdb_code, chain_id):
        path = CONFIGS.SCRATCH_OUT + pdb_code + chain_id + CONFIGS.DOT_ACC
        return self.get_one_hot(path, CONFIGS.ACC2)

scratch = Scratch()
print(scratch.get_secondary_structure("6y2d", "A"))
print(scratch.get_solvency_accessibility("6y2d", "A"))