#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Time    : 2022/1/10 12:22
# @Author  : zhangbc0315@outlook.com
# @File    : rcts_config.py
# @Software: PyCharm

from config.config_base import ConfigBase


class RCTSConfig(ConfigBase):

    is_retro = False
    name = "centralized"
    col_num = "c_num"
    col_code = "rct_code"

    def __init__(self, config_json, rxn_data_tsv_fp: str, rct_tsv_fp: str, rxn_code_rt_code_tsv_fp: str, evaluate_year: int,
                 cat_code_and_ohi_fp: str, sol_code_and_ohi_fp: str):
        self.rxn_data_tsv_fp = rxn_data_tsv_fp
        self.rxn_template_tsv_fp = rct_tsv_fp
        self.rxn_code_rt_code_tsv_fp = rxn_code_rt_code_tsv_fp
        self.evaluate_year = evaluate_year
        self.cat_code_and_ohi_fp = cat_code_and_ohi_fp
        self.sol_code_and_ohi_fp = sol_code_and_ohi_fp

        super(RCTSConfig, self).__init__(config_json)

        self.filter_tids_fp = self.get_file_path_by_key("filter_tids_file_name")
        self.train_rxn_code_and_rt_code_fp = self.get_file_path_by_key("train_rids_file_name")
        self.test_rxn_code_and_rt_code_fp = self.get_file_path_by_key("test_rids_file_name")
        self.train_temp_dp = self.get_file_path_by_key("train_temp_dir_name")
        self.test_temp_dp = self.get_file_path_by_key("test_temp_dir_name")
        self.model_dp = self.get_file_path_by_key("model_dir_name", create_dir=True)

        self.min_num_covered_rxns_by_rxn_template = config_json.get("min_num_covered_rxns_by_rxn_template")
        self.device = config_json.get("device")
        self.batch_size = config_json.get("batch_size")
        self.epoch_num = config_json.get("epoch_num")
        self.lr_start = config_json.get("lr_start")
        self.lr_end = config_json.get("lr_end")


if __name__ == "__main__":
    pass
