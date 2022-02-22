#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Time    : 2021/8/25 14:48
# @Author  : zhangbc0315@outlook.com
# @File    : prepare_batch_data.py
# @Software: PyCharm

import sys
import os
from random import shuffle

import torch
from torch_geometric.data import Batch

from config import RDPath


class PrepareBatchData:

    # _right_single_dp = RDPath.RXN_RIGHT_SINGLE_DATA_DP
    # _wrong_single_dp = RDPath.RXN_WRONG_SINGLE_DATA_DP
    # _train_batch_dp = RDPath.RXN_TRAIN_BATCH_GNN_DATA_DP
    # _test_batch_dp = RDPath.RXN_TEST_BATCH_GNN_DATA_DP

    _right_single_dp = "C:\\Users\\zhang\\Documents\\Data\\RxnMolPredictor\\RxnDeterminer\\rxn_right_single_data"
    _wrong_single_dp = "C:\\Users\\zhang\\Documents\\Data\\RxnMolPredictor\\RxnDeterminer\\rxn_wrong_single_data"
    _train_batch_dp = "C:\\Users\\zhang\\Documents\\Data\\RxnMolPredictor\\RxnDeterminer\\rxn_train_batch_gnn_data"
    _test_batch_dp = "C:\\Users\\zhang\\Documents\\Data\\RxnMolPredictor\\RxnDeterminer\\rxn_test_batch_gnn_data"

    @classmethod
    def random_fns(cls):
        right_fns = list(os.listdir(cls._right_single_dp))
        wrong_fns = list(os.listdir(cls._wrong_single_dp))
        fns = right_fns + wrong_fns
        shuffle(fns)

        length = len(fns)
        train_fns = fns[: int(0.9*length)]
        test_fns = fns[int(0.9*length):]
        print(f"总文件: {len(fns)}, 训练文件: {len(train_fns)}, 测试文件: {len(test_fns)}")
        return train_fns, test_fns

    @classmethod
    def save_fns(cls):
        print(f'saving {RDPath.RXN_TRAIN_SINGLE_FNS_FP}, {RDPath.RXN_TEST_SINGLE_FNS_FP}')
        train_fns, test_fns = cls.random_fns()
        with open(RDPath.RXN_TRAIN_SINGLE_FNS_FP, 'w', encoding='utf-8')as f:
            f.write('\n'.join(train_fns))
        with open(RDPath.RXN_TEST_SINGLE_FNS_FP, 'w', encoding='utf-8')as f:
            f.write('\n'.join(test_fns))

    @classmethod
    def load_fns(cls):
        print(f'loading {RDPath.RXN_TRAIN_SINGLE_FNS_FP}, {RDPath.RXN_TEST_SINGLE_FNS_FP}')
        with open(RDPath.RXN_TRAIN_SINGLE_FNS_FP, 'r', encoding='utf-8')as f:
            train_fns = f.read().split('\n')
        with open(RDPath.RXN_TEST_SINGLE_FNS_FP, 'r', encoding='utf-8')as f:
            test_fns = f.read().split('\n')
        print(f'loaded')
        return train_fns, test_fns

    @classmethod
    def get_fns(cls):
        if not os.path.exists(RDPath.RXN_TEST_SINGLE_FNS_FP) or not os.path.exists(RDPath.RXN_TRAIN_SINGLE_FNS_FP):
            cls.save_fns()
        return cls.load_fns()

    @classmethod
    def save_batch_data_list(cls, r_data_list, p_data_list, batch_id, batch_dp):
        r_fp = os.path.join(batch_dp, f"r_{batch_id}.pt")
        p_fp = os.path.join(batch_dp, f"p_{batch_id}.pt")
        r_batch = Batch.from_data_list(r_data_list)
        p_batch = Batch.from_data_list(p_data_list)
        torch.save(r_batch, r_fp)
        torch.save(p_batch, p_fp)
        print(f"save batch {batch_id}")

    @classmethod
    def get_batch_data(cls, fns: [str], batch_dp: str, pid):
        r_data_list = []
        p_data_list = []
        batch_id = 0
        for n, fn in enumerate(fns):
            if (n // 1000) % 5 != pid:
                continue
            if len(r_data_list) == 1000:
                cls.save_batch_data_list(r_data_list, p_data_list, batch_id*5+pid, batch_dp)
                batch_id += 1
                r_data_list = []
                p_data_list = []
            y = 1 if fn.startswith('r') else 0
            dp = cls._right_single_dp if fn.startswith('r') else cls._wrong_single_dp
            fp = os.path.join(dp, fn)
            single_data = torch.load(fp)
            r_data = single_data['rs']
            p_data = single_data['p']

            cat_len = single_data['cat'].size()[0]
            sol_len = single_data['sol'].size()[0]
            r_data.g = torch.reshape(torch.cat([single_data['cat'], single_data['sol']]), [1, cat_len + sol_len])
            r_data.y = y

            r_data_list.append(r_data)
            p_data_list.append(p_data)
        if len(r_data_list) != 0:
            cls.save_batch_data_list(r_data_list, p_data_list, batch_id*5 + pid, batch_dp)

    @classmethod
    def process(cls, pid: int):
        train_fns, test_fns = cls.get_fns()
        cls.get_batch_data(train_fns, cls._train_batch_dp, pid)
        cls.get_batch_data(test_fns, cls._test_batch_dp, pid)


if __name__ == "__main__":
    i = int(sys.argv[1])
    PrepareBatchData.process(i)
