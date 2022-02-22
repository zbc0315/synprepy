#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Time    : 2022/2/18 15:52
# @Author  : zhangbc0315@outlook.com
# @File    : retro_rcts_config.py
# @Software: PyCharm

from config.rcts_config import RCTSConfig


class RetroRCTSConfig(RCTSConfig):

    is_retro = True
    name = "centralized"
    col_num = "c_num"
    col_code = "rct_code"

    def __init__(self, config_json, rxn_data_tsv_fp: str, rct_tsv_fp: str, rxn_code_rt_code_tsv_fp: str, evaluate_year: int,
                 cat_code_and_ohi_fp: str, sol_code_and_ohi_fp: str):
        super(RetroRCTSConfig, self).__init__(config_json, rxn_data_tsv_fp, rct_tsv_fp, rxn_code_rt_code_tsv_fp,
                                              evaluate_year, cat_code_and_ohi_fp, sol_code_and_ohi_fp)


if __name__ == "__main__":
    pass
