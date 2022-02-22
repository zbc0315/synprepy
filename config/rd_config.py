#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Time    : 2022/2/22 11:34
# @Author  : zhangbc0315@outlook.com
# @File    : rd_config.py
# @Software: PyCharm

from config.config_ml import ConfigML


class RDConfig(ConfigML):

    def __init__(self, config_json, rxn_data_fp: str, rct_fp: str, evaluate_year: int):
        self.rxn_data_fp = rxn_data_fp
        self.evaluate_year = evaluate_year
        self.rct_fp = rct_fp
        super(RDConfig, self).__init__(config_json)

        self.low_yield_rxns_tsv_fp = self.get_file_path_by_key("low_yield_rxns_tsv_file_name")


if __name__ == "__main__":
    pass
