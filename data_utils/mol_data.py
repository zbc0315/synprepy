#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Time    : 2022/2/11 14:56
# @Author  : zhangbc0315@outlook.com
# @File    : mol_data.py
# @Software: PyCharm

from data_utils.mol import Mol

import pandas as pd


class MolData:

    def __init__(self, mol_fp: str):
        self._mol_df: pd.DataFrame = pd.read_csv(mol_fp, sep='\t', encoding='utf-8')

    def get_mol_by_mol_code(self, mol_code: str):
        mol_df = self._mol_df.query(f"mol_code=='{mol_code}'")
        if len(mol_df) == 0:
            return None
        mol_df = mol_df.reset_index()
        return Mol(mol_code, mol_df["smiles"][0], mol_df["inchi"][0])


if __name__ == "__main__":
    pass
