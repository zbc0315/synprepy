#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Time    : 2022/1/4 11:01
# @Author  : zhangbc0315@outlook.com
# @File    : add_ys.py
# @Software: PyCharm
import os

import torch
import pandas as pd

from data_utils.tid_utils import TidUtils


class AddYs:

    @classmethod
    def _process(cls, dp, tsv_fp):
        df = pd.read_csv(tsv_fp, sep='\t', encoding='utf-8')
        for idx, row in df.iterrows():
            if idx % 1000 == 0:
                print(idx)
            fp = os.path.join(dp, f"{idx}.pt")
            data = torch.load(fp)
            if data.key != row.product_smi_with_map:
                raise ValueError(f"smiles is wrong: {idx}")
            encoded_tid = TidUtils.encode_tids([row.tid])
            data.y = torch.tensor([encoded_tid], dtype=torch.float)
            torch.save(data, fp)

    @classmethod
    def process(cls):
        cls._process("C:\\Users\\zhang\\Documents\\Data\\RxnMolPredictor\\BaseTemp\\consider_bond_energy\\train_single_data\\",
                     "C:\\Users\\zhang\\OneDrive\\Projects\\RXN-PREDICTOR\\Data\\DataBase\\rxndb_ml_data_train_gnn.tsv",)
        cls._process(
            "C:\\Users\\zhang\\Documents\\Data\\RxnMolPredictor\\BaseTemp\\consider_bond_energy\\test_single_data\\",
            "C:\\Users\\zhang\\OneDrive\\Projects\\RXN-PREDICTOR\\Data\\DataBase\\rxndb_ml_data_test_gnn.tsv", )


if __name__ == "__main__":
    AddYs.process()
