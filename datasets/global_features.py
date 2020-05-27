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