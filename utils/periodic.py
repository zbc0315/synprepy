#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Time    : 2022/1/8 21:44
# @Author  : zhangbc0315@outlook.com
# @File    : periodic.py
# @Software: PyCharm

import os

import pandas as pd


class Periodic:

    _kernels = [86, 54, 36, 18, 10, 2]

    def __init__(self):
        periodic_fp = "periodic.tsv"
        if not os.path.exists(periodic_fp):
            periodic_fp = "utils/periodic.tsv"
        if not os.path.exists(periodic_fp):
            periodic_fp = "../utils/periodic.tsv"
        if not os.path.exists(periodic_fp):
            raise FileNotFoundError(f"找不到文件 periodic.tsv，请检查 utils/periodic.tsv 是否存在。")
        self._df = pd.read_csv(periodic_fp, sep='\t', encoding='utf-8')

    def get_row_and_col(self, atomic_num: int):
        row = 1
        col = 1
        for i, k in enumerate(self._kernels):
            if atomic_num > k:
                row = 7 - i
                col = atomic_num - k
                break
        return row, col

    def get_row(self, atomic_num: int) -> int:

        return self._df.loc[atomic_num-1, 'row']

    def get_col(self, atomic_num: int) -> int:
        return self._df.loc[atomic_num-1, 'column']


if __name__ == "__main__":
    p = Periodic()
    print(p.get_col(10))
    print(p.get_row(10))
