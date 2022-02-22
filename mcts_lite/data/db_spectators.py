#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Time    : 2021/9/3 11:06
# @Author  : zhangbc0315@outlook.com
# @File    : db_spectators.py
# @Software: PyCharm

import pandas as pd

from mcts_lite.mcts_path import MctsPath


class DbSpectators:

    cats_df = pd.read_csv(MctsPath.CATS_COUNT_FILTER_FP, sep='\t', encoding='utf-8')
    sols_df = pd.read_csv(MctsPath.SOLS_COUNT_FILTER_FP, sep='\t', encoding='utf-8')
    cats_id_code_df = pd.read_csv(MctsPath.CATS_IDX_TO_WEBMID_FP, sep='\t', encoding='utf-8')
    sols_id_code_df = pd.read_csv(MctsPath.SOLS_IDX_TO_WEBMID_FP, sep='\t', encoding='utf-8')

    @classmethod
    def cats_idxes_to_codes(cls, cat_idxes):
        res = []
        for idx in cat_idxes:
            code = cls.cat_id_to_code(idx)
            if code is not None:
                res.append(code)
        return res

    @classmethod
    def cat_id_to_code(cls, idx):
        if idx < 0 or idx >= len(cls.cats_id_code_df):
            return None
        return cls.cats_id_code_df.loc[idx, 'mol_code']

    @classmethod
    def sols_idxes_to_codes(cls, sol_idxes):
        res = []
        for idx in sol_idxes:
            code = cls.sol_id_to_code(idx)
            if code is not None:
                res.append(code)
        return res

    @classmethod
    def sol_id_to_code(cls, idx):
        if idx < 0 or idx >= len(cls.sols_id_code_df):
            return None
        return cls.sols_id_code_df.loc[idx, 'mol_code']


if __name__ == "__main__":
    pass
