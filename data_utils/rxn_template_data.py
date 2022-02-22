#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Time    : 2022/2/15 23:00
# @Author  : zhangbc0315@outlook.com
# @File    : rxn_template_data.py
# @Software: PyCharm

from typing import Optional
import logging

import pandas as pd

from data_utils.rxn_template import RxnTemplate


class RxnTemplateData:

    def __init__(self, tsv_fp: str):
        self._rt_df = pd.read_csv(tsv_fp, sep='\t', encoding='utf-8')
        logging.info(f"load {tsv_fp}")

    def init_filter(self, min_num_covered_by_rxn: int):
        logging.info(f"before filter, there is {len(self._rt_df)} rt")
        self._rt_df = self._rt_df.query(f"num_of_covered_rxn>={min_num_covered_by_rxn}")
        logging.info(f"after filter by {min_num_covered_by_rxn}, there is {len(self._rt_df)} rt")

    def __len__(self) -> int:
        return len(self._rt_df)

    def get_all_template(self):
        for i, row in self._rt_df.iterrows():
            rt = RxnTemplate(row.rt_code,
                             row.rxn_template,
                             row.num_of_covered_rxn)
            yield rt

    def get_rt_smiles_by_rt_code(self, rt_code: str) -> Optional[int]:
        query = self._rt_df.query(f"rt_code=='{rt_code}'")
        if len(query) == 0:
            return None
        query = query.reset_index()
        return query.rxn_template[0]

    def get_rt_smiles_by_ohi(self, ohi: int) -> str:
        return self._rt_df["rxn_template"][ohi]

    def get_one_hot_idx_by_rt_code(self, rt_code: str) -> Optional[int]:
        query = self._rt_df.query(f"rt_code=='{rt_code}'")
        if len(query) == 0:
            return None
        return query.index[0]

    def get_template_by_rt_smiles(self, smiles: str) -> Optional[RxnTemplate]:
        res = self._rt_df.query(f"rxn_template=='{smiles}'")
        if len(res) == 0:
            return None
        res = res.reset_index()
        return RxnTemplate(res.rt_code[0], res.rxn_template[0], res.num_of_covered_rxn[0])


if __name__ == "__main__":
    pass
