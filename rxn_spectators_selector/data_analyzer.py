# -*- coding: utf-8 -*-
# @Time     : 2021/4/27 9:47
# @Author   : zbc
# @Email    : zhangbc0315@outlook.com
# @FileName : data_analyzer.py
# @Software : PyCharm

import json

import numpy as np
import pandas as pd

from config import SSPath


class DataAnalyzer:

    _cats_count = {"smiles": [], "times": []}
    _sols_count = {"smiles": [], "times": []}

    @classmethod
    def _add_smiles(cls, dic: {}, smiles: str):
        if smiles not in dic["smiles"]:
            dic["smiles"].append(smiles)
            dic["times"].append(1)
        else:
            idx = dic["smiles"].index(smiles)
            dic["times"][idx] += 1

    @classmethod
    def _count_spectators(cls, role: str, df: pd.DataFrame):
        if role == "cats":
            dic = cls._cats_count
        else:
            dic = cls._sols_count
        for i in df.index:
            if i % 100000 == 0:
                print(f"{role} - {i}")
            spectators = df.loc[i, role]
            if spectators is None or not isinstance(spectators, str):
                continue
            spectators_list = json.loads(spectators)
            if len(spectators_list) == 0:
                smiles = 0
            else:
                smiles = spectators_list[0]
            cls._add_smiles(dic, smiles)

    @classmethod
    def process(cls):
        df = pd.read_csv(SSPath.RXN_WITH_SPECTATORS_FP, sep="\t", encoding="utf-8")
        cls._count_spectators("cats", df)
        cls._count_spectators("sols", df)
        cats_df = pd.DataFrame(cls._cats_count)
        sols_df = pd.DataFrame(cls._sols_count)
        cats_df.sort_values("times", inplace=True)
        sols_df.sort_values("times", inplace=True)
        cats_df["ratio"] = cats_df["times"] / cats_df["times"].sum()
        sols_df["ratio"] = sols_df["times"] / sols_df["times"].sum()
        cats_df.to_csv(SSPath.CATS_COUNT_FP, sep='\t', index=False, encoding="utf-8")
        sols_df.to_csv(SSPath.SOLS_COUNT_FP, sep='\t', index=False, encoding="utf-8")


if __name__ == '__main__':
    DataAnalyzer.process()
