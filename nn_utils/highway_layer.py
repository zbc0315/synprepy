#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Time    : 2021/8/10 12:47
# @Author  : zhangbc0315@outlook.com
# @File    : highway_layer.py
# @Software: PyCharm

import torch
import torch.nn as nn
import torch.nn.functional as F


class Highway(nn.Module):

    def __init__(self, input_size: int):
        super(Highway, self).__init__()
        self.relu = nn.ReLU()
        self.sigmoid = nn.Sigmoid()
        self.proj = nn.Linear(input_size, input_size)
        self.transform = nn.Linear(input_size, input_size)
        self.transform.bias.data.fill_(-2.0)

    def forward(self, data):
        proj_result = self.relu(self.proj(data))
        proj_gate = self.sigmoid(self.transform(data))
        gated = (proj_gate * proj_result) + ((1 - proj_gate) * data)
        return gated


if __name__ == "__main__":
    pass
