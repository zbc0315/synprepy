#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Time    : 2022/2/12 15:40
# @Author  : zhangbc0315@outlook.com
# @File    : yp_config.py
# @Software: PyCharm

from config.config_base import ConfigBase


class YPConfig(ConfigBase):
    
    def __init__(self, config_json,
                 mol_data_fp: str,
                 rxn_data_fp: str,
                 cat_code_to_ohi_fp: str,
                 sol_code_to_ohi_fp: str,
                 evaluate_year: int):
        self.mol_data_fp = mol_data_fp
        self.rxn_data_fp = rxn_data_fp
        self.cat_code_to_ohi_fp = cat_code_to_ohi_fp
        self.sol_code_to_ohi_fp = sol_code_to_ohi_fp
        self.evaluate_year = evaluate_year
        super(YPConfig, self).__init__(config_json)
        self.train_rxn_codes_fp = self.get_file_path_by_key("train_rxn_codes_file_name")
        self.test_rxn_codes_fp = self.get_file_path_by_key("test_rxn_codes_file_name")
        self.train_temp_dp = self.get_file_path_by_key("train_temp_dir_name", create_dir=True)
        self.test_temp_dp = self.get_file_path_by_key("test_temp_dir_name", create_dir=True)
        self.model_dp = self.get_file_path_by_key("model_dir_name", create_dir=True)
        self.device = config_json.get("device")
        self.batch_size = config_json.get("batch_size")
        self.epoch_num = config_json.get("epoch_num")
        self.lr_start = config_json.get("lr_start")
        self.lr_end = config_json.get("lr_end")


if __name__ == "__main__":
    pass
