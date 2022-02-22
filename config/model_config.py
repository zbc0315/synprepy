#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Time    : 2022/2/22 11:02
# @Author  : zhangbc0315@outlook.com
# @File    : model_config.py
# @Software: PyCharm

import json


class ModelConfig:

    def __init__(self, config_fp: str):
        self._config_fp = config_fp
        with open(config_fp, 'r', encoding='utf-8')as f:
            self._config_json = json.load(f)
        self._num_input = self._config_json.get("nn_num_input")

    def _dump_config(self):
        with open(self._config_fp, 'w', encoding='utf-8') as f:
            self._config_json.dump(f)

    # region ===== num_node_features =====
    @property
    def gnn_num_node_features(self):
        return self._config_json.get("gnn_num_node_features")

    @gnn_num_node_features.setter
    def gnn_num_node_features(self, value):
        self._config_json["gnn_num_node_features"] = value
        self._dump_config()
    # endregion

    # region ===== gnn_num_edge_features =====
    @property
    def gnn_num_edge_features(self):
        return self._config_json.get("gnn_num_edge_features")

    @gnn_num_edge_features.setter
    def gnn_num_edge_features(self, value):
        self._config_json["gnn_num_edge_features"] = value
        self._dump_config()
    # endregion

    # region ===== gnn_num_global_features =====
    @property
    def gnn_num_global_features(self):
        return self._config_json.get("gnn_num_global_features")

    @gnn_num_global_features.setter
    def gnn_num_global_features(self, value):
        self._config_json["gnn_num_global_features"] = value
        self._dump_config()
    # endregion

    # region ===== num_out_channels =====
    @property
    def num_out_channels(self):
        return self._config_json.get("num_out_channels")

    @num_out_channels.setter
    def num_out_channels(self, value):
        self._config_json["num_out_channels"] = value
        self._dump_config()
    # endregion

    # region ===== nn_num_features =====
    @property
    def nn_num_features(self):
        return self._config_json.get("nn_num_features")

    @nn_num_features.setter
    def nn_num_features(self, value):
        self._config_json["nn_num_features"] = value
        self._dump_config()
    # endregion


if __name__ == "__main__":
    pass
