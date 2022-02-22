# -*- coding: utf-8 -*-
# @Time     : 2021/5/26 16:26
# @Author   : zbc
# @Email    : zhangbc0315@outlook.com
# @FileName : prepare_right_and_high_yield_rxn.py
# @Software : PyCharm

import pandas as pd

from rmp_database.old_chem_db import OldChemDb
from config import RDPath


class PrepareRightAndHighYieldRxn:

    @classmethod
    def _load_right_rxns(cls):
        df = pd.read_csv(RDPath.RIGHT_RXN_FP, sep='\t', encoding="utf-8")
        return df

    @classmethod
    def _load_yield_rxn(cls):
        odb = OldChemDb()
        for rxn_data in odb.get_data_iter("public.yield_rxn", ["rid", "product_inchi", "percent_yield"], None):
            yield rxn_data

    @classmethod
    def _get_right_and_yield_rxn(cls):
        right_rxn_df = cls._load_right_rxns()
        res_dict = {"rid": [], "cat_smiles": [], "sol_smiles": [], "product_inchi": [], "percent_yield": []}
        for i, yield_rxn_data in enumerate(cls._load_yield_rxn()):
            if i % 10000 == 0:
                print(i)
            rxn_spec = right_rxn_df[right_rxn_df["rid"] == yield_rxn_data["rid"]]
            if rxn_spec["rid"].count() == 0:
                continue
            rxn_spec.reset_index(inplace=True)
            res_dict["rid"].append(yield_rxn_data["rid"])
            res_dict["cat_smiles"].append(rxn_spec["cat"][0])
            res_dict["sol_smiles"].append(rxn_spec["sol"][0])
            res_dict["product_inchi"].append(yield_rxn_data["product_inchi"])
            res_dict["percent_yield"].append(yield_rxn_data["percent_yield"])
        res_df = pd.DataFrame(res_dict)
        res_df.to_csv(RDPath.RIGHT_AND_YIELD_RXN_FP, sep='\t', encoding='utf-8', index=False)

    @classmethod
    def _get_right_and_high_yield_rxn(cls):
        right_and_yield_rxn_df = pd.read_csv(RDPath.RIGHT_AND_YIELD_RXN_FP, sep='\t')
        right_and_high_yield_rxn_df = right_and_yield_rxn_df[right_and_yield_rxn_df["percent_yield"] >= 85]
        right_and_high_yield_rxn_df.to_csv(RDPath.RIGHT_AND_HIGH_YIELD_RXN_FP, sep='\t', encoding="utf-8", index=False)

    @classmethod
    def process(cls):
        cls._get_right_and_yield_rxn()
        cls._get_right_and_high_yield_rxn()


if __name__ == '__main__':
    PrepareRightAndHighYieldRxn.process()
