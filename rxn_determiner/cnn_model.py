# -*- coding: utf-8 -*-
# @Time     : 2021/6/8 15:30
# @Author   : zbc
# @Email    : zhangbc0315@outlook.com
# @FileName : cnn_model.py
# @Software : PyCharm

from torch import nn
from torch.nn import Module


class SampleCnn(Module):

    def __init__(self, input_dim, output_dim):
        super().__init__()
        self.dense_layer = nn.Linear(input_dim, 126)


if __name__ == '__main__':
    pass
