# -*- coding: utf-8 -*-
# @Time     : 2021/3/10 14:33
# @Author   : zbc
# @Email    : zhangbc0315@outlook.com
# @FileName : global_config.py
# @Software : PyCharm

import os
import json


class GlobalConfig:

    _config = None

    @classmethod
    def _load_config(cls):
        if cls._config is not None:
            return
        with open("setting.json", 'r', encoding='utf-8')as f:
            cls._config = json.load(f)

    @classmethod
    def _init_proj(cls):
        cls._load_config()
        if not os.path.exists(cls._config['proj_path']):
            os.mkdir(cls._config['proj_path'])

    @classmethod
    def _init_base_template_selector(cls):
        pass


if __name__ == '__main__':
    pass
