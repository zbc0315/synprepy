#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Time    : 2022/2/11 18:18
# @Author  : zhangbc0315@outlook.com
# @File    : rxn_data_test.py
# @Software: PyCharm

import unittest

from config.config import Config
from data_utils.rxn_data import MolData, RxnData, RxnDataType


class RxnDataTest(unittest.TestCase):

    def test_query_rxn(self):
        config = Config("../../config.json")
        mol_data = MolData(config.mol_data_tsv_file_path)
        rxn_data = RxnData(config.rxn_data_tsv_file_path, mol_data, RxnDataType.TRAIN_TEST, config.evaluate_year)
        rxn = rxn_data.get_rxn_by_rxn_code("R2141291")
        self.assertEqual(rxn.rxn_mols_code, "1439613>7452,14>1439614")

    def test_query_rxn_2(self):
        config = Config("../../config.json")
        mol_data = MolData(config.mol_data_tsv_file_path)
        rxn_data = RxnData(config.rxn_data_tsv_file_path, mol_data, RxnDataType.TRAIN_TEST, config.evaluate_year)
        rxn = rxn_data.get_rxn_by_rxn_code("R285534")
        self.assertIsNone(rxn)


if __name__ == "__main__":
    unittest.main()
