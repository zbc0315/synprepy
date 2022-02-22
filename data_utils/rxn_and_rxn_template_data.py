#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Time    : 2022/1/10 12:41
# @Author  : zhangbc0315@outlook.com
# @File    : rxn_and_rxn_template_data.py
# @Software: PyCharm

from typing import Optional

import pandas as pd


class RxnAndRxnTemplateData:

    _rt_type_to_keys = {"centralized": ("c_tid", "c_num"),
                        "extended": ("e_tid", "e_num")}

    def __init__(self, rxn_fp: str, rt_fp: str, rid_tid_fp: str, template_type: str):
        self._rxn_df = pd.read_csv(rxn_fp, sep='\t', encoding='utf-8')
        self._rt_df = pd.read_csv(rt_fp, sep='\t', encoding='utf-8')
        self._rid_tid_df = pd.read_csv(rid_tid_fp, sep='\t', encoding='utf-8')

        self._tid_col, self._num_covered_col = self._rt_type_to_keys[template_type]

    def get_rids_and_tids_by_covered(self, min_num_covered: int) -> ([str], [int]):
        res_df = self._rid_tid_df.query(f"{self._num_covered_col}>={min_num_covered}")
        return res_df['rxn_code'], res_df[self._tid_col]

    def get_rxn_smiles_by_rid(self, rid: int) -> Optional[str]:
        res_df = self._rxn_df.query(f"rid=={rid}")
        if len(res_df) == 0:
            return None
        res_df = res_df.reset_index()
        return res_df.loc[0, 'rxn_smi']

    def get_tid_by_rid(self, rid: int) -> Optional[int]:
        res_df = self._rid_tid_df.query(f"rid=={rid}")
        if len(res_df) == 0:
            return None
        res_df = res_df.reset_index()
        return res_df.loc[0, self._tid_col]

    def get_template_smiles_by_tid(self, tid: int) -> Optional[str]:
        res_df = self._rt_df.query(f"tid=={tid}")
        if len(res_df) == 0:
            return None
        res_df = res_df.reset_index()
        return res_df.loc[0, "rxn_template"]


if __name__ == "__main__":
    pass
