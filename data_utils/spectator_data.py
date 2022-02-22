#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Time    : 2022/2/12 13:21
# @Author  : zhangbc0315@outlook.com
# @File    : spectator_data.py
# @Software: PyCharm

from enum import Enum
import logging

import pandas as pd


class SpectatorType(Enum):
    COMMON = 0,
    RARE = 1,
    NONE = 2


class SpectatorData:

    def __init__(self, cat_code_to_one_hot_fp: str):
        self._cat_code_to_one_hot_idx_df = pd.read_csv(cat_code_to_one_hot_fp, sep='\t', encoding='utf-8')
        logging.info(f"load: {cat_code_to_one_hot_fp}")

    @property
    def num_common(self) -> int:
        return len(self._cat_code_to_one_hot_idx_df)

    @property
    def len_one_hot(self):
        return self.num_common + 2

    @property
    def rare_one_hot_idx(self) -> int:
        return self.num_common

    @property
    def none_one_hot_idx(self) -> int:
        return self.num_common + 1

    def get_type_by_one_hot_idx(self, one_hot_idx: int):
        if 0 <= one_hot_idx < self.num_common:
            return SpectatorType.COMMON
        elif one_hot_idx == self.rare_one_hot_idx:
            return SpectatorType.RARE
        elif one_hot_idx == self.none_one_hot_idx:
            return SpectatorType.NONE
        else:
            raise ValueError(f"expect 0 <= one_hot_idx < {self.len_one_hot}, but get {one_hot_idx}")

    def get_one_hot_idx_by_mol_code(self, mol_code: str) -> int:
        if mol_code is None:
            return self.none_one_hot_idx
        spec_df = self._cat_code_to_one_hot_idx_df.query(f"mol_code=='{mol_code}'")
        if len(spec_df) == 0:
            return self.rare_one_hot_idx
        spec_df = spec_df.reset_index()
        return spec_df["one_hot_idx"][0]

    def get_one_hot_idx_by_mol_codes(self, mol_codes: [str]) -> int:
        if len(mol_codes) == 0:
            return self.none_one_hot_idx
        elif len(mol_codes) > 1:
            return self.rare_one_hot_idx
        else:
            return self.get_one_hot_idx_by_mol_code(mol_codes[0])

    def get_one_hot_vec_by_mol_codes(self, mol_codes: [str]) -> [int]:
        one_hot_idx = self.get_one_hot_idx_by_mol_codes(mol_codes)
        ohv = [0]*self.len_one_hot
        ohv[one_hot_idx] = 1
        return ohv

    def get_one_hot_vec_by_mol_code(self, mol_code: str) -> [int]:
        one_hot_idx = self.get_one_hot_idx_by_mol_code(mol_code)
        ohv = [0]*self.len_one_hot
        ohv[one_hot_idx] = 1
        return ohv


if __name__ == "__main__":
    pass
