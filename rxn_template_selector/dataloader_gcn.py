#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Time    : 2021/8/24 14:22
# @Author  : zhangbc0315@outlook.com
# @File    : dataloader_gcn.py
# @Software: PyCharm

import os

import torch

from config import BTPath, CTPath


class DataloaderGcn:
    batch_data_dps = {"base_temp": {"train": BTPath.CBE_BATCH_DATA_TRAIN_DP, "test": BTPath.CBE_BATCH_DATA_TEST_DP},
                      "comp_temp": {"train": CTPath.GCN_TRAIN_BATCH_DATA_DP, "test": CTPath.GCN_TEST_BATCH_DATA_DP}}
    batch_data_dps = {"base_temp": {"test": "C:\\Users\\zhang\\Documents\\Data\\RxnPredictor\\GNN\\BaseTemplate_AtomicNum\\train_batch_1000_data",
                                    "train": "C:\\Users\\zhang\\Documents\\Data\\RxnPredictor\\GNN\\BaseTemplate_AtomicNum\\test_batch_1000_data"},
                      "comp_temp": {'train': 'C:\\Users\\zhang\\Documents\\Data\\RxnMolPredictor\\CompTemp\\gcn\\train_batch_1000_data',
                                    'test': 'C:\\Users\\zhang\\Documents\\Data\\RxnMolPredictor\\CompTemp\\gcn\\test_batch_1000_data'}}
    # batch_data_dps = {"base_temp": {"train": "E:\\C\\Documents\\Data\\RxnPredictor\\GNN\\BaseTemplate_AtomicNum\\train_processed_batch_1000",
    #                                 "test": "E:\\C\\Documents\\Data\\RxnPredictor\\GNN\\BaseTemplate_AtomicNum\\test_processed_batch_1000"}}
    single_data_dps = {"base_temp": BTPath.CBE_SINGLE_DATA_TEST_DP,
                       'comp_temp': CTPath.GCN_SINGLE_DATA_DP}

    def __init__(self, proj: str, data_type: str):
        self.proj = proj
        self.data_dp = self.batch_data_dps[proj][data_type]
        print(f"load data from: {self.data_dp}")
        self.data_fns = list(os.listdir(self.data_dp))

    def __getitem__(self, item):
        return torch.load(os.path.join(self.data_dp, self.data_fns[item]))

    def __len__(self):
        return len(self.data_fns)

    def get_data_param(self):
        single_data_fp = self.single_data_dps[self.proj]
        fp = os.path.join(single_data_fp, '0.pt')
        data = torch.load(fp)
        res = {"num node features": data.x.size()[1],
               "num edge features": data.edge_attr.size()[1],
               "num out channels": data.y.size().numel()}
        return res


if __name__ == "__main__":
    dg = DataloaderGcn('base_temp', 'train')
    print(dg.get_data_param())
