#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Time    : 2022/1/5 12:54
# @Author  : zhangbc0315@outlook.com
# @File    : rxn_handler_test.py
# @Software: PyCharm

import unittest

from rdkit.Chem import AllChem

from chem_utils.chem_handler import RxnHandler


class RxnHandlerTest(unittest.TestCase):

    def test_remove_unmapped_mols_in_rxn1(self):
        rxn_smarts = "[CH3:1]CCC.CCCCC>>[CH3:1]CCCCCCC"
        rxn = AllChem.ReactionFromSmarts(rxn_smarts)
        rxn = RxnHandler.remove_unmapped_mols_in_rxn(rxn)
        new_rxn_smarts = AllChem.ReactionToSmiles(rxn)
        self.assertEqual(new_rxn_smarts, "CCC[CH3:1]>>CCCCCCC[CH3:1]")


if __name__ == "__main__":
    unittest.main()
