#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Time    : 2021/9/13 11:06
# @Author  : zhangbc0315@outlook.com
# @File    : prepare_data.py
# @Software: PyCharm
import pandas as pd

from rmp_database.old_chem_db import OldChemDb
from mct.mcts_path import MctsPath


class PrepareData:

    _odb = OldChemDb()

    @classmethod
    def _get_products_data(cls):
        odb = OldChemDb()
        for rxn_data in odb.get_data_iter('rxn_for_test.clean_duplicate_rxn', ['rid', 'rxn_code'], f'p_num=1'):
            yield rxn_data['rid'], rxn_data['rxn_code'].split('>')[-1]

    @classmethod
    def _search_mol_inchi(cls, mid):
        mol_data = cls._odb.search_one_data('rxn_for_test.mols', ['inchi'], f'mid={mid}')
        return mol_data['inchi']

    @classmethod
    def process(cls):
        res_dict = {'rid': [], 'inchi': []}
        for n, (rid, p_mid) in enumerate(cls._get_products_data()):
            p_inchi = cls._search_mol_inchi(p_mid)
            res_dict['rid'].append(rid)
            res_dict['inchi'].append(p_inchi)
            if n % 10000 == 0:
                print(n)
        df = pd.DataFrame(res_dict)
        df.to_csv(MctsPath.TARGET_PRODUCTS_FP, sep='\t', encoding='utf-8', index=False)


if __name__ == "__main__":
    PrepareData.process()
