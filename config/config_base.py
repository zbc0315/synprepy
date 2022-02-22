#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Time    : 2022/1/10 14:11
# @Author  : zhangbc0315@outlook.com
# @File    : config_base.py
# @Software: PyCharm
import abc
import os
import logging


class ConfigBase:

    def __init__(self, config_json: {}):
        self._config_json = config_json
        self._root = self._config_json.get("root")
        if not os.path.join(self._root):
            logging.info(f"Create dir: {self._root}")
            os.mkdir(self._root)

    def get_file_path_by_key(self, key: str, check_exists: bool = False, create_dir: bool = False) -> str:
        fn = self._config_json.get(key)
        fp = os.path.join(self._root, fn)
        if check_exists and not os.path.exists(fp):
            logging.warning(f"Can't find {fp} for {key}")
        if create_dir and not os.path.exists(fp):
            os.mkdir(fp)
        return fp


if __name__ == "__main__":
    pass
