#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Time    : 2022/2/11 18:09
# @Author  : zhangbc0315@outlook.com
# @File    : mol_data_test.py
# @Software: PyCharm

import unittest

from config.config import Config
from data_utils.mol_data import Mol, MolData


class MolDataTest(unittest.TestCase):

    def test_query_mol(self):
        config = Config("../../config.json")
        mol_data = MolData(config.mol_data_tsv_file_path)
        mol = mol_data.get_mol_by_mol_code("M5")
        self.assertEqual(mol.inchi, "InChI=1S/CH2O3.Ag/c2-1(3)4;/h(H2,2,3,4);/q;+2/p-2")


if __name__ == "__main__":
    unittest.main()
