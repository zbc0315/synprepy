#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Time    : 2022/2/19 10:47
# @Author  : zhangbc0315@outlook.com
# @File    : nn_base_model.py
# @Software: PyCharm

import torch.nn as nn
import torch.nn.functional as F


class Highway(nn.Module):

    def __init__(self, input_size):
        super(Highway, self).__init__()
        self.elu = nn.ELU()
        self.proj = nn.Linear(input_size, input_size)
        self.transform = nn.Linear(input_size, input_size)
        self.transform.bias.data.fill_(-2.0)

    def forward(self, in_data):
        proj_result = self.elu(self.proj(in_data))
        proj_gate = self.elu(self.transform(in_data))
        gated = (proj_gate * proj_result) + ((1 - proj_gate) * in_data)
        return gated


if __name__ == "__main__":
    pass
