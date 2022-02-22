#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Time    : 2021/8/16 16:19
# @Author  : zhangbc0315@outlook.com
# @File    : dataloader_highway.py
# @Software: PyCharm

import os

import torch

from config import BTPath
from config import CTPath


class DataloaderHighway:

    # data_dps = {"base_temp": {"train": BTPath.TRAIN_BATCH_DATA_DP, "test": BTPath.TEST_BATCH_DATA_DP},
    #             "comp_temp": {"train": CTPath.HIGHWAY_TRAIN_BATCH_DATA_DP, "test": CTPath.HIGHWAY_TEST_BATCH_DATA_DP}}
    data_dps = {"base_temp": {"train": "C:\\Users\\zhang\\Documents\\Data\\RxnMolPredictor\\BaseTemp\\train_batch_data_1000",
                              "test": "C:\\Users\\zhang\\Documents\\Data\\RxnMolPredictor\\BaseTemp\\test_batch_data_1000"},
                "comp_temp": {"train": "C:\\Users\\zhang\\Documents\\Data\\RxnMolPredictor\\CompTemp\\highway_gcn\\train_batch_1000_data",
                              "test": "C:\\Users\\zhang\\Documents\\Data\\RxnMolPredictor\\CompTemp\\highway_gcn\\test_batch_1000_data"}}

    def __init__(self, proj: str = "base_temp", data_type: str = "train"):
        self.data_dp = self.data_dps[proj][data_type]
        self.data_fns = list(os.listdir(self.data_dp))

    def __len__(self):
        return len(self.data_fns)

    def __getitem__(self, item: int):
        return torch.load(os.path.join(self.data_dp, self.data_fns[item]))


if __name__ == "__main__":
    pass
