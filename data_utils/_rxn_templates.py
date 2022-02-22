#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Time    : 2022/1/7 14:07
# @Author  : zhangbc0315@outlook.com
# @File    : _rxn_templates.py
# @Software: PyCharm

import pandas as pd


class RxnTemplates:

    def __init__(self, data_fp: str):
        self._df = pd.read_csv(data_fp, sep='\t', encoding='utf-8')
        self._rxn_tot = None

    def query_by_tid(self, tid: int):
        res = self._df.query(f"tid == {tid}")
        if len(res) == 0:
            return None
        res = res.reset_index()
        if len(res) > 1:
            # TODO 修改为logger
            print(f"WARN: 存在多个tid={tid}")
        return res.iloc[0]

    def query_by_temp(self, temp: str):
        res = self._df.query(f"rxn_template=='{temp}'")
        if len(res) == 0:
            return None
        res = res.reset_index()
        if len(res) > 1:
            raise ValueError(f"Can't contain more than one same template, but get: {len(res)}")
        return res.iloc[0]

    def get_max_covered(self):
        return self._df.num_of_covered_rxn.max()

    def calc_coverage(self, min_num_covered: int):
        if self._rxn_tot is None:
            self._rxn_tot = self._df.num_of_covered_rxn.sum()
        filter_df = self._df.query(f"num_of_covered_rxn>={min_num_covered}")
        return filter_df.num_of_covered_rxn.sum() / self._rxn_tot


if __name__ == "__main__":
    pass
