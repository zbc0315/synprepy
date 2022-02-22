#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Time    : 2021/9/2 14:59
# @Author  : zhangbc0315@outlook.com
# @File    : db_rxn_template.py
# @Software: PyCharm
import math

import pandas as pd

from mcts_lite.mcts_path import MctsPath


class DbRxnTemplate:

    print('loading rxn template db')
    basic_df = pd.read_csv(MctsPath.BASIC_RXN_TEMP_FP, sep='\t', encoding='utf-8')
    comp_df = pd.read_csv(MctsPath.COMP_RXN_TEMP_FP, sep='\t', encoding='utf-8')
    print('loaded rxn template db')

    @classmethod
    def search_basic_smates_by_tid(cls, tid: int) -> (str, float):
        if tid == 35800:
            tid = 97000
        smates, rid_num = cls._search_smates_by_tid(tid, cls.basic_df)
        cover_score = rid_num / 85000
        # cover_score = math.sqrt(rid_num) / 400
        cover_score = 0
        # if rid_num < 90:
        #     return None, None
        if tid in [97400, 65900, 99800, 1722000]:
            return None, None
        else:
            return smates, cover_score

    @classmethod
    def search_comp_smates_by_tid(cls, tid: int):
        smates, rid_num = cls._search_smates_by_tid(tid, cls.comp_df)
        return smates

    @classmethod
    def _search_smates_by_tid(cls, tid: int, df: pd.DataFrame):
        res_df = df.query(f'tid=={tid}')
        if len(res_df.index) == 0:
            return None
        return res_df.smates[res_df.index[0]], res_df.rid_num[res_df.index[0]]


if __name__ == "__main__":
    print(DbRxnTemplate.search_comp_smates_by_tid(1066310))
