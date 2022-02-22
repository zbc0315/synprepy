#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Time    : 2022/2/21 13:45
# @Author  : zhangbc0315@outlook.com
# @File    : cs_ss_config.py
# @Software: PyCharm

from config.config_ml import ConfigML


class CSSSConfig(ConfigML):

    def __init__(self, config_json, rxn_data_tsv_fp: str, evaluate_year: int, spe_code_to_ohi_fp: str,
                 rxn_code_and_spe_code_fp: str):
        self.rxn_data_tsv_fp = rxn_data_tsv_fp
        self.evaluate_year = evaluate_year
        self.spectator_code_to_ohi_fp = spe_code_to_ohi_fp
        self.rxn_code_and_spe_code_fp = rxn_code_and_spe_code_fp

        super(CSSSConfig, self).__init__(config_json)


if __name__ == "__main__":
    pass
