#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Time    : 2022/2/15 22:58
# @Author  : zhangbc0315@outlook.com
# @File    : rxn_template.py
# @Software: PyCharm


class RxnTemplate:

    def __init__(self, rt_code: str, rt_smiles: str, count: int):
        self._rt_code = rt_code
        self._rt_smiles = rt_smiles
        self._count = count

    # region ===== rt_code =====
    @property
    def rt_code(self):
        return self._rt_code

    @rt_code.setter
    def rt_code(self, value):
        self._rt_code = value
    # endregion

    # region ===== rt_smiles =====
    @property
    def rt_smiles(self):
        return self._rt_smiles

    @rt_smiles.setter
    def rt_smiles(self, value):
        self._rt_smiles = value
    # endregion

    # region ===== count =====
    @property
    def count(self):
        return self._count

    @count.setter
    def count(self, value):
        self._count = value
    # endregion


if __name__ == "__main__":
    pass
