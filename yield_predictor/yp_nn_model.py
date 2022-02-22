#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Time    : 2022/2/20 10:38
# @Author  : zhangbc0315@outlook.com
# @File    : yp_nn_model.py
# @Software: PyCharm
import torch
from torch.nn import Module, Linear
from torch.nn import functional as F

from utils.nn_base_model import Highway


class YPNNModel(Module):

    def __init__(self, input_size: int, spectators_size: int):
        super(YPNNModel, self).__init__()
        self.lin_r1 = Linear(input_size + spectators_size, 256)
        self.lin_r2 = Linear(256, 256)
        self.lin_p1 = Linear(input_size, 256)
        self.lin_p2 = Linear(256, 256)

        self.hr1 = Highway(256)
        self.hr2 = Highway(256)
        self.hr3 = Highway(256)
        self.hr4 = Highway(256)
        self.hr5 = Highway(256)

        self.lin = Linear(512, 1)

    def forward(self, r_data, p_data):
        r = F.elu(self.lin_r1(r_data))
        r = F.elu(self.lin_r2(r))
        r = F.elu(self.hr1(r))
        r = F.dropout(r, 0.2)
        r = F.elu(self.hr2(r))
        r = F.dropout(r, 0.2)
        r = F.elu(self.hr3(r))
        r = F.dropout(r, 0.2)
        r = F.elu(self.hr4(r))
        r = F.dropout(r, 0.2)
        r = F.elu(self.hr5(r))

        p = F.elu(self.lin_p1(p_data))
        p = F.elu(self.lin_p2(p))

        # return F.cosine_similarity(r, p)*100
        return F.sigmoid(self.lin(torch.cat([r, p], dim=1))) * 100


if __name__ == "__main__":
    pass
