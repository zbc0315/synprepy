#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Time    : 2021/8/25 17:16
# @Author  : zhangbc0315@outlook.com
# @File    : dataloader_gcn.py
# @Software: PyCharm

import os

import torch

from config import RDPath


class DataloaderGcn:

    # dps = {'train': RDPath.RXN_TRAIN_BATCH_GNN_DATA_DP,
    #        'test': RDPath.RXN_TEST_BATCH_GNN_DATA_DP}

    dps = {'train': 'C:\\Users\\zhang\\Documents\\Data\\RxnMolPredictor\\RxnDeterminer\\rxn_train_batch_gnn_data',
           'test': 'C:\\Users\\zhang\\Documents\\Data\\RxnMolPredictor\\RxnDeterminer\\rxn_test_batch_gnn_data'}

    def __init__(self, train_test: str):
        self.dp = self.dps[train_test]
        self.fns = list(os.listdir(self.dp))

    def __getitem__(self, item):
        r_fp = os.path.join(self.dp, f'r_{item}.pt')
        p_fp = os.path.join(self.dp, f'p_{item}.pt')

        return torch.load(r_fp), torch.load(p_fp)

    def __len__(self):
        return len(self.fns) // 2

    @classmethod
    def get_data_param(cls):
        single_data_fp = os.path.join(RDPath.RXN_RIGHT_SINGLE_DATA_DP, 'r_0.pt')
        data = torch.load(single_data_fp)
        res = {'num node features': data['rs'].x.size()[1],
               'num edge features': data['rs'].edge_attr.size()[1],
               'num global features': data['cat'].size()[0] + data['sol'].size()[0]}
        return res


if __name__ == "__main__":
    print(DataloaderGcn.get_data_param())
