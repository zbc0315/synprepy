#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Time    : 2021/9/22 11:36
# @Author  : zhangbc0315@outlook.com
# @File    : config.py
# @Software: PyCharm

import os
import json
import socket


class Config:

    _config_json = None
    _host_name = socket.gethostname()

    @classmethod
    def consider_spectators(cls):
        return cls.config_json()['consider_spectators']

    @classmethod
    def device(cls):
        return cls._get_config_by_host_name('device')

    @classmethod
    def data_dp(cls):
        return cls._get_config_by_host_name('data_dp')

    # =======================================

    @classmethod
    def _get_config_by_host_name(cls, config_name):
        configs = cls.config_json()[config_name]
        if cls._host_name in configs.keys():
            return configs[cls._host_name]
        return configs['default']

    # =======================================

    @classmethod
    def config_json(cls):
        if cls._config_json is None:
            config_path = cls._get_config_file_path()
            cls._load_config_fp(config_path)
        return cls._config_json

    # =====================================

    @classmethod
    def _load_config_fp(cls, config_fp: str):
        with open(config_fp, 'r', encoding='utf-8')as f:
            cls._config_json = json.load(f)

    @classmethod
    def _get_config_file_path(cls):
        config_fp = 'config.json'
        for i in range(3):
            if os.path.exists(config_fp):
                print(f'config fp: {config_fp}')
                return config_fp
            config_fp = os.path.join('..', config_fp)


if __name__ == "__main__":
    print(Config.config_json())
