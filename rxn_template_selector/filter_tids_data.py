#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Time    : 2022/1/10 13:37
# @Author  : zhangbc0315@outlook.com
# @File    : filter_tids_data.py
# @Software: PyCharm

import os
import logging

import pandas as pd


class FilterTidsData:

    def __init__(self, filter_tids_fp: str, min_num_covered: int = None, template_fp: str = None):
        self._filter_tids_fp = filter_tids_fp
        if not os.path.exists(filter_tids_fp):
            if min_num_covered is None and template_fp is None:
                raise AttributeError(f"Expect min_num_covered and template_fp are not None "
                                     f"while {filter_tids_fp} not exists, "
                                     f"but got min_num_covered or template_fp is None")
            self._get_filter_tids(template_fp, min_num_covered)
        self._filter_tids_df = pd.read_csv(filter_tids_fp, sep='\t', encoding='utf-8')

    @property
    def num_tids(self) -> int:
        return len(self._filter_tids_df)

    def get_one_hot_idx_by_tid(self, tid: int) -> int:
        res_df = self._filter_tids_df.query(f"tid=={tid}")
        if len(res_df) == 0:
            raise ValueError(f"Can't find tid: {tid} in {self._filter_tids_fp}")
        res_df = res_df.reset_index()
        return res_df.loc[0, "one_hot_idx"]

    def _get_filter_tids(self, template_fp: str, min_num_covered):
        logging.info(f"Filter tids to {template_fp} with min num covered {min_num_covered}")
        template_df = pd.read_csv(template_fp, sep='\t', encoding='utf-8')
        filter_df = template_df.query(f"num_of_covered_rxn>={min_num_covered}")
        filter_df['one_hot_idx'] = range(len(filter_df))
        filter_df.to_csv(self._filter_tids_fp, sep='\t', encoding='utf-8', index=False)


if __name__ == "__main__":
    pass
