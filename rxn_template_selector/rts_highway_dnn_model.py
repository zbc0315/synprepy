#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Time    : 2021/8/13 16:51
# @Author  : zhangbc0315@outlook.com
# @File    : rts_highway_dnn_model.py
# @Software: PyCharm

from typing import Tuple
import os

import torch
import torch.nn as nn
from torch.nn import Linear, Dropout
import torch.nn.functional as F

from nn_utils.highway_layer import Highway
from config import BTPath


class HighwayDnnModel(nn.Module):

    def __init__(self, dim_input, dim_out):
        super(HighwayDnnModel, self).__init__()
        # self.lin1 = nn.Linear(dim_input, 512)
        # self.lin2 = nn.Linear(512, dim_out)
        self.elu = nn.ELU()
        self.linear1 = Linear(dim_input, 512)
        self.dropout1 = Dropout(0.3)
        self.highways = []
        self.dropouts = []
        self.h1 = Highway(512)
        self.h2 = Highway(512)
        self.h3 = Highway(512)
        self.h4 = Highway(512)
        self.h5 = Highway(512)
        self.d1 = Dropout(0.1)
        self.d2 = Dropout(0.1)
        self.d3 = Dropout(0.1)
        self.d4 = Dropout(0.1)
        self.d5 = Dropout(0.1)
        self.linear2 = Linear(512, dim_out)

    def forward(self, data):
        x = self.elu(self.linear1(data))
        x = self.dropout1(x)
        x = self.d1(self.elu(self.h1(x)))
        x = self.d2(self.elu(self.h2(x)))
        x = self.d3(self.elu(self.h3(x)))
        x = self.d4(self.elu(self.h4(x)))
        x = self.d5(self.elu(self.h5(x)))
        x = self.linear2(x)
        return x
        # y = F.relu(self.lin1(data))
        # return self.lin2(y)

    @classmethod
    def load_model(cls, dim_in, dim_out, device: str, model_dp: str):
        fns = list(os.listdir(model_dp))
        if len(fns) == 0:
            print(f"没有已保存的模型")
            return cls(dim_in, dim_out).to(device), 0, None
        fns = sorted(fns, key=lambda x: float(x.split('_')[2]))
        best_fp = os.path.join(model_dp, fns[-1])

        es = [int(fn.split('_')[0]) for fn in fns]
        max_e = max(es)

        print(f"加载模型：{best_fp}")
        model = cls(dim_in, dim_out).to(device)
        # model = torch.load(best_fp).to(device)
        best_state = torch.load(best_fp)
        model.load_state_dict(best_state['model'])
        return model, max_e, best_state['opt']


if __name__ == "__main__":
    pass
