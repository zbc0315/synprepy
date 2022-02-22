#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Time    : 2021/9/15 13:45
# @Author  : zhangbc0315@outlook.com
# @File    : get_spectators_to_webmid.py
# @Software: PyCharm

import pandas as pd
from rdkit.Chem import AllChem

from mcts_lite.mcts_path import MctsPath
from synpre_backend.webdb.mol_table import MolTable


class GetSpectatorsToWebmid:

    @classmethod
    def _get_id_and_smiles(cls, source_fp):
        df = pd.read_csv(source_fp, sep='\t', encoding='utf-8', index_col=0)
        for i, row in df.iterrows():
            yield i, row.smiles

    @classmethod
    def _process(cls, source_fp, target_fp):
        res = {'mol_idx': [], 'mol_code': []}
        for idx, smiles in cls._get_id_and_smiles(source_fp):
            mol = AllChem.MolFromSmiles(smiles)
            mol_code = MolTable.save_mol(mol)
            res['mol_idx'].append(idx)
            res['mol_code'].append(mol_code)
        res_df = pd.DataFrame(res)
        print(f'save to file: {target_fp}')
        res_df.to_csv(target_fp, sep='\t', encoding='utf-8', index=False)

    @classmethod
    def process(cls):
        cls._process(MctsPath.CATS_COUNT_FILTER_FP, MctsPath.CATS_IDX_TO_WEBMID_FP)
        cls._process(MctsPath.SOLS_COUNT_FILTER_FP, MctsPath.SOLS_IDX_TO_WEBMID_FP)


if __name__ == "__main__":
    GetSpectatorsToWebmid.process()
