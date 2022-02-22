#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Time    : 2021/10/22 14:30
# @Author  : zhangbc0315@outlook.com
# @File    : get_cat_sol_name.py
# @Software: PyCharm

import pandas as pd

from synpre_backend.webdb.mol_table import MolTable


CAT_IDX_TO_MOL_CODE_FP = 'C:\\Users\\zhang\\Documents\\Data\\RxnMolPredictor\\Mcts\\CatPredictor\\cats_idx_to_webmid.tsv'
SOL_IDX_TO_MOL_CODE_FP = 'C:\\Users\\zhang\\Documents\\Data\\RxnMolPredictor\\Mcts\\SolPredictor\\sols_idx_to_webmid.tsv'


class GetCatSolName:

    cat_idx_to_mol_code_df = pd.read_csv(CAT_IDX_TO_MOL_CODE_FP, sep='\t', encoding='utf-8')
    sol_idx_to_mol_code_df = pd.read_csv(SOL_IDX_TO_MOL_CODE_FP, sep='\t', encoding='utf-8')

    @classmethod
    def get_mol_code_from_cat_idx(cls, cat_idx: int):
        line = cls.cat_idx_to_mol_code_df.query(f'mol_idx=={cat_idx}')
        return line.mol_code[line.index[0]]

    @classmethod
    def get_mol_code_from_sol_idx(cls, sol_idx: int):
        line = cls.sol_idx_to_mol_code_df.query(f'mol_idx=={sol_idx}')
        return line.mol_code[line.index[0]]

    @classmethod
    def process_cat_idx(cls, cat_idx):
        print('-' * 100)
        if cat_idx == -1:
            print('CAT: NONE')
        elif cat_idx == 325:
            print('CAT: OTHERS')
        else:
            mol_code = cls.get_mol_code_from_cat_idx(cat_idx)
            mol = MolTable.get_mol_by_mol_code(mol_code)
            print(f'CAT: {cat_idx}')
            print(f'{cat_idx}: {mol.molName}')
            print(mol.smiles)
            # print(mol)

    @classmethod
    def process_sol_idx(cls, sol_idx):
        print('-'*100)
        if sol_idx == -1:
            print('SOL: NONE')
        elif sol_idx == 321:
            print('SOL: OTHERS')
        else:
            mol_code = cls.get_mol_code_from_sol_idx(sol_idx)
            mol = MolTable.get_mol_by_mol_code(mol_code)
            print(f'SOL: {sol_idx}')
            print(f'{sol_idx}: {mol.molName}')
            print(mol.smiles)
            # print(mol)

    @classmethod
    def process_cats_idxes(cls, cats_idxes):
        for cat_idx in cats_idxes.split('.'):
            cat_idx = int(cat_idx)
            cls.process_cat_idx(cat_idx)

    @classmethod
    def process_sols_idxes(cls, sols_idxes):
        for sol_idx in sols_idxes.split('.'):
            sol_idx = int(sol_idx)
            cls.process_sol_idx(sol_idx)


if __name__ == "__main__":
    print('='*200)
    GetCatSolName.process_cats_idxes('323.319.317.320.305.321.285.271.-1.325')
    print('='*200)
    GetCatSolName.process_sols_idxes('312.315.318.302.282.320.314.313.-1.321')

