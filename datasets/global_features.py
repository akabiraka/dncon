import sys
sys.path.append("../dncon")
import numpy as np
import pandas as pd
import torch
import torch.nn.functional as F

class GlobalFeatures(object):
    def __init__(self):
        super(GlobalFeatures, self).__init__()

    def encode_protein_length(self, l):
        """
        encode protein length (<75: 1000, 75–150: 0100, 150–225: 0010, >225: 0001)
        """
        values = 0
        if(l<=75):
            values = np.full(l, 0)
        elif l>75 and l<=150:
            values = np.full(l, 1)
        elif l>150 and l<=225:
            values = np.full(l, 2)
        else:
            values = np.full(l, 3)
        one_hot_tensor = F.one_hot(torch.tensor(values), num_classes=4)
        return pd.DataFrame(one_hot_tensor.numpy())

    def encode_window_position(self, pos, l):
        """
        pos: current window position
        l: protein length
        """
        relative_pos = pos/l
        feature = np.zeros(l)
        feature[:] = relative_pos
        return pd.DataFrame(feature)

    def __get_numeric_value(self, pos):
        aa_separation_intervals = [1, 13, 19, 27, 39, 51, 63, 75, 87, 99, 111, 121]
        for i, value in enumerate(aa_separation_intervals):
            if i == 11:
                return i
            elif pos >= aa_separation_intervals[i] and pos < aa_separation_intervals[i+1]:
                return i

    def encode_residue_separation(self, aa_pos1, aa_pos2):
        """
        aa_pos1, aa_pos2: the positions of two amino-acid in protein sequence
        """
        encoded_pos1 = self.__get_numeric_value(aa_pos1)
        encoded_pos2 = self.__get_numeric_value(aa_pos2)
        encoded_pos1 = F.one_hot(torch.tensor(encoded_pos1), num_classes=12)
        encoded_pos2 = F.one_hot(torch.tensor(encoded_pos2), num_classes=12)
        return encoded_pos1 | encoded_pos2

    def encode_residue_separation_df(self, aa_pos1, aa_pos2):
        one_hot_tensor = self.encode_residue_separation(aa_pos1, aa_pos2)
        return pd.DataFrame(one_hot_tensor.numpy())

gf = GlobalFeatures()
# gf.encode_window_position(5, 155)
print(gf.encode_residue_separation_df(10, 28))