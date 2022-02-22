#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Time    : 2022/2/11 14:58
# @Author  : zhangbc0315@outlook.com
# @File    : mol.py
# @Software: PyCharm


class Mol:

    def __init__(self, mol_code: str, smiles: str, inchi: str):
        self._mol_code: str = mol_code
        self._smiles: str = smiles
        self._inchi: str = inchi

    # region ===== mol_code =====
    @property
    def mol_code(self):
        return self._mol_code

    @mol_code.setter
    def mol_code(self, value):
        self._mol_code = value
    # endregion

    # region ===== smiles =====
    @property
    def smiles(self):
        return self._smiles

    @smiles.setter
    def smiles(self, value):
        self._smiles = value
    # endregion

    # region ===== inchi =====
    @property
    def inchi(self):
        return self._inchi

    @inchi.setter
    def inchi(self, value):
        self._inchi = value
    # endregion


if __name__ == "__main__":
    pass
