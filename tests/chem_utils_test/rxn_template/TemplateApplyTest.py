#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Time    : 2022/2/18 21:53
# @Author  : zhangbc0315@outlook.com
# @File    : TemplateApplyTest.py
# @Software: PyCharm

import unittest

from rdkit.Chem import AllChem

from chem_utils.rxn_template.template_apply import TemplateApply


class TemplateApplyTest(unittest.TestCase):

    def test_template_apply1(self):
        temp_smi = "(Cl[c:1].[NH2:2])>>[c:1][NH:2]"
        # temp_smi = "(Cl[C:2].[NH2:1])>>[NH:1][C:2]"
        p_smi = "CSc1nc(Cl)cc(Nc2ccc(C(=O)OC(C)(C)C)cc2)n1"
        res = TemplateApply.retro_synthesis_by_smarts_to_smiles(temp_smi, p_smi)
        self.assertEqual(res, [['CC(C)(C)OC(=O)c1ccc(N)cc1.CSc1nc(Cl)cc(Cl)n1'], ['CC(C)(C)OC(=O)c1ccc(Cl)cc1.CSc1nc(N)cc(Cl)n1']])
        # print(res)


if __name__ == "__main__":
    unittest.main()
