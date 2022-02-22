#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Time    : 2022/1/6 15:53
# @Author  : zhangbc0315@outlook.com
# @File    : config.py
# @Software: PyCharm

import os
import json
import logging

from config.yp_config import YPConfig
from config.retro_rcts_config import RetroRCTSConfig
from config.rcts_config import RCTSConfig
from config.rets_config import RETSConfig
from config.config_base import ConfigBase
from config.cs_ss_config import CSSSConfig
from config.rd_config import RDConfig


class Config(ConfigBase):

    def __init__(self, config_json_fp):
        with open(config_json_fp, 'r', encoding='utf-8')as f:
            _config_json = json.load(f)
        super(Config, self).__init__(_config_json)

        self.evaluate_year = _config_json["evaluate_year"]
        self.rxn_data_tsv_file_path = self.get_file_path_by_key("rxn_data_tsv_file_name", True)
        self.mol_data_tsv_file_path = self.get_file_path_by_key("mol_data_tsv_file_name", True)
        self.rxn_code_with_rxn_template_tsv_file_path = self.get_file_path_by_key("rxn_code_with_rxn_template_tsv_file_name")
        self.rxn_centralized_template_tsv_file_path = self.get_file_path_by_key("rxn_centralized_template_tsv_file_name")
        self.rxn_extended_template_tsv_file_path = self.get_file_path_by_key("rxn_extended_template_tsv_file_name")
        self.rxn_code_with_rt_code_tsv_file_path = self.get_file_path_by_key("rxn_code_with_rt_code_tsv_file_name")
        self.min_num_covered_rxns_by_cat = _config_json["min_num_covered_rxns_by_cat"]
        self.min_num_covered_rxns_by_sol = _config_json["min_num_covered_rxns_by_sol"]
        self.cat_code_and_one_hot_idx_tsv_file_path = self.get_file_path_by_key("cat_code_and_one_hot_idx_tsv_file_name")
        self.sol_code_and_one_hot_idx_tsv_file_path = self.get_file_path_by_key("sol_code_and_one_hot_idx_tsv_file_name")
        self.rxn_code_and_cat_code_fp = self.get_file_path_by_key("rxn_code_and_cat_code_tsv_file_name")
        self.rxn_code_and_sol_code_fp = self.get_file_path_by_key("rxn_code_and_sol_code_tsv_file_name")

        self.rcts_config = RCTSConfig(_config_json.get("rxn_centralized_template_selector_config"),
                                      self.rxn_data_tsv_file_path,
                                      self.rxn_centralized_template_tsv_file_path,
                                      self.rxn_code_with_rt_code_tsv_file_path,
                                      self.evaluate_year,
                                      self.cat_code_and_one_hot_idx_tsv_file_path,
                                      self.sol_code_and_one_hot_idx_tsv_file_path)
        self.retro_rcts_config = RetroRCTSConfig(_config_json.get("retro_rxn_centralized_template_selector_config"),
                                                 self.rxn_data_tsv_file_path,
                                                 self.rxn_centralized_template_tsv_file_path,
                                                 self.rxn_code_with_rt_code_tsv_file_path,
                                                 self.evaluate_year,
                                                 self.cat_code_and_one_hot_idx_tsv_file_path,
                                                 self.sol_code_and_one_hot_idx_tsv_file_path)
        self.rets_config = RETSConfig(_config_json.get("rxn_extended_template_selector_config"),
                                      self.rxn_data_tsv_file_path,
                                      self.rxn_extended_template_tsv_file_path,
                                      self.rxn_code_with_rt_code_tsv_file_path,
                                      self.evaluate_year,
                                      self.cat_code_and_one_hot_idx_tsv_file_path,
                                      self.sol_code_and_one_hot_idx_tsv_file_path)
        self.rd_config = RDConfig(_config_json.get("rxn_determiner_config"),
                                  self.rxn_data_tsv_file_path,
                                  self.rxn_centralized_template_tsv_file_path,
                                  self.evaluate_year)
        self.yp_config = YPConfig(_config_json.get("yield_predictor_config"),
                                  self.mol_data_tsv_file_path,
                                  self.rxn_data_tsv_file_path,
                                  self.cat_code_and_one_hot_idx_tsv_file_path,
                                  self.sol_code_and_one_hot_idx_tsv_file_path,
                                  self.evaluate_year)
        self.cs_config = CSSSConfig(_config_json.get("cat_selector_config"),
                                    self.rxn_data_tsv_file_path,
                                    self.evaluate_year,
                                    self.cat_code_and_one_hot_idx_tsv_file_path,
                                    self.rxn_code_and_cat_code_fp)
        self.ss_config = CSSSConfig(_config_json.get("sol_selector_config"),
                                    self.rxn_data_tsv_file_path,
                                    self.evaluate_year,
                                    self.sol_code_and_one_hot_idx_tsv_file_path,
                                    self.rxn_code_and_sol_code_fp)


if __name__ == "__main__":
    pass
