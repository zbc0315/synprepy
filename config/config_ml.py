#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Time    : 2022/2/21 13:41
# @Author  : zhangbc0315@outlook.com
# @File    : config_ml.py
# @Software: PyCharm

from config.config_base import ConfigBase


class ConfigML(ConfigBase):

    def __init__(self, config_json):
        super(ConfigML, self).__init__(config_json)

        self.train_rxn_code_fp = self.get_file_path_by_key("train_rxn_codes_file_name")
        self.test_rxn_code_fp = self.get_file_path_by_key("test_rxn_codes_file_name")
        self.train_temp_dp = self.get_file_path_by_key("train_temp_dir_name")
        self.test_temp_dp = self.get_file_path_by_key("test_temp_dir_name")
        self.model_dp = self.get_file_path_by_key("model_dir_name", create_dir=True)

        self.device = config_json.get("device")
        self.batch_size = config_json.get("batch_size")
        self.epoch_num = config_json.get("epoch_num")
        self.lr_start = config_json.get("lr_start")
        self.lr_end = config_json.get("lr_end")


if __name__ == "__main__":
    pass
