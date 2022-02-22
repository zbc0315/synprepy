#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Time    : 2021/9/26 21:33
# @Author  : zhangbc0315@outlook.com
# @File    : count_tids.py
# @Software: PyCharm

import os

import torch

from config import BTPath


class CountTids:

    train_batch_dp = "C:\\Users\\zhang\\Documents\\Data\\BaseTemp\\consider_bond_energy\\train_batch_1000_data"

    @classmethod
    def parse_y(cls, ys):
        return ys.max(dim=1)[1]

    @classmethod
    def count(cls):
        tids = set()
        for n, batch_fn in enumerate(os.listdir(cls.train_batch_dp)):
            batch_fp = os.path.join(cls.train_batch_dp, batch_fn)
            batch = torch.load(batch_fp)
            new_tids = cls.parse_y(batch.y)
            tids = tids.union(set(new_tids.tolist()))
            if n % 100 == 0:
                print(f'{n} - {len(tids)}')
        tids = sorted(tids, key=lambda x: x)
        tids = [str(tid) for tid in tids]
        with open(BTPath.CBE_TIDS_FP, 'w', encoding='utf-8')as f:
            f.write('\n'.join(tids))


if __name__ == "__main__":
    CountTids.count()
