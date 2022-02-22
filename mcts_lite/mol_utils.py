#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Time    : 2021/9/13 13:50
# @Author  : zhangbc0315@outlook.com
# @File    : mol_utils.py
# @Software: PyCharm

import pandas as pd

from mcts_lite.data.db_cac import DbCac


class MolUtils:

    mol_df = None

    @classmethod
    def save(cls, fp: str):
        if cls.mol_df is not None:
            cls.mol_df.to_csv(fp, sep='\t', encoding='utf-8')

    @classmethod
    def init(cls):
        cls.mol_df = pd.DataFrame(columns=['inchi', 'is_cac'])

    @classmethod
    def add_or_query_inchis(cls, inchis):
        for inchi in inchis:
            yield cls.add_or_query_inchi(inchi)

    @classmethod
    def add_or_query_inchi(cls, inchi):
        idx = cls._query_id_by_inchi(inchi)
        if idx is None:
            idx = cls._add_inchi(inchi)
        return idx

    @classmethod
    def query_by_mid(cls, mid):
        return cls.mol_df.loc[mid, ]

    @classmethod
    def _add_inchi(cls, inchi):
        is_cac = DbCac.is_cac(inchi=inchi)
        cls.mol_df = cls.mol_df.append({'inchi': inchi, 'is_cac': is_cac}, ignore_index=True)
        return len(cls.mol_df) - 1

    @classmethod
    def _query_id_by_inchi(cls, inchi):
        if len(cls.mol_df) == 0:
            return None
        res = cls.mol_df.query(f'inchi=="{inchi}"')
        if len(res) == 0:
            return None
        idx = res.index[0]
        return idx


if __name__ == "__main__":
    DbCac.init('123')
    MolUtils.init()
    id1 = MolUtils.add_or_query_inchi('123')
    id2 = MolUtils.add_or_query_inchi('123')
    id3 = MolUtils.add_or_query_inchi('234')
    id4 = MolUtils.add_or_query_inchi('234')
    r0 = MolUtils.query_by_mid(0)
    r1 = MolUtils.query_by_mid(1)
    print('1')
