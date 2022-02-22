#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Time    : 2022/2/20 11:37
# @Author  : zhangbc0315@outlook.com
# @File    : dataset.py
# @Software: PyCharm

from torch_geometric import data as gdata


class Dataset(gdata.InMemoryDataset):

    def __init__(self, train_or_test: str, rt_config):
        super(Dataset, self).__init__()


if __name__ == "__main__":
    pass
