# -*- coding: utf-8 -*-
# @Time     : 2021/5/31 16:41
# @Author   : zbc
# @Email    : zhangbc0315@outlook.com
# @FileName : prepare_right_rxn.py
# @Software : PyCharm

import json

import pandas as pd

from config import RDPath
from rmp_database.old_chem_db import OldChemDb


class PrepareRightRxn:

    @classmethod
    def _is_right_spectators(cls, spectators_str: str):
        if not isinstance(spectators_str, str):
            return False, ""
        spectators = json.loads(spectators_str)
        spectators = list(set(spectators))
        if len(spectators) == 0:
            return True, "0"
        elif len(spectators) == 1:
            return True, spectators[0]
        else:
            return False, ""

    @classmethod
    def _get_rxn_with_spectators(cls):
        rxn_df = pd.read_csv(RDPath.RXN_WITH_SPECTATORS_FP, sep='\t', encoding="utf-8")
        for i, row in rxn_df.iterrows():
            cats = row["cats"]
            sols = row["sols"]
            cat_is_right, cat_smiles = cls._is_right_spectators(cats)
            sol_is_right, sol_smiles = cls._is_right_spectators(sols)
            if not (cat_is_right and sol_is_right):
                continue
            yield row["rid"], cat_smiles, sol_smiles

    @classmethod
    def process(cls):
        cats = set()
        sols = set()
        res_dic = {"rid": [], "cat": [], "sol": []}
        for i, (rid, cat_smi, sol_smi) in enumerate(cls._get_rxn_with_spectators()):
            if i % 10000 == 0:
                print(i)
            res_dic["rid"].append(rid)
            res_dic["cat"].append(cat_smi)
            res_dic["sol"].append(sol_smi)
        res_df = pd.DataFrame(res_dic)
        res_df.to_csv(RDPath.RIGHT_RXN_FP, sep='\t', encoding="utf-8", index=False)


if __name__ == '__main__':
    PrepareRightRxn.process()
