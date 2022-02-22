#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Time    : 2021/9/14 10:03
# @Author  : zhangbc0315@outlook.com
# @File    : crkit_tester.py
# @Software: PyCharm

from rdkit.Chem import AllChem

from crkitpy.inchiparse.inchi_parser import InchiParser
from crkitpy.smatesparse.smates_parser import SmatesParser
from crkitpy.rxnhandler.rxn_template_apply import RxnTemplateApply
from crkitpy.inchiparse.inchi_writer import InchiWriter
from crkitpy.smilesparse.smiles_writer import SmilesWriter


class CrkitTester:

    @classmethod
    def get_product_and_temp(cls):
        product_inchi = 'InChI=1S/C23H22FN5O/c1-15-12-28(13-16-5-3-2-4-6-16)14-19(15)21-26-23(30)20-11-25-22(29(20)27-21)17-7-9-18(24)10-8-17/h2-11,15,19H,12-14H2,1H3,(H,26,27,30)'
        # print(AllChem.MolToSmiles(AllChem.MolFromInchi(product_inchi)))
        smates = '[NH2:1][cH2:5][cH2:4][NH:3][CH:6]=[O:2]>>[nH:1]1[cH:5][cH:4][nH:3][cH:6]1'
        # product_inchi = 'InChI=1S/C17H23N3O2/c1-17(2,3)22-16(21)20(4)11-7-8-13-12(10-11)15-14(19-13)6-5-9-18-15/h5-6,9,11,19H,7-8,10H2,1-4H3'
        # smates = '[CH3:1][CH:7]=[O:4].[I:2][cH2:6][cH2:5][NH2:3]>>[cH:1]1[cH:6][cH:5][nH:3][cH:7]1'
        # product_inchi = 'InChI=1S/C15H15F3N4O3/c1-7-5-10(14(24)25-2)20-12-9(6-19-22(7)12)13(23)21-11(8-3-4-8)15(16,17)18/h5-6,8,11H,3-4H2,1-2H3,(H,21,23)'
        # smates = '[NH3:1].[O:2]=[CH2:3]>>[NH:1]=[CH2:3]'
        # product_inchi = 'InChI=1S/C18H18N4O2/c1-21-13-9-5-3-7-11(13)19-17(21)15(23)16(24)18-20-12-8-4-6-10-14(12)22(18)2/h3-10,15-16,23-24H,1-2H3'
        # smates = '[Br:2][CH2:5][CH:6]=[O:3].[NH2:1][CH:7]=[NH:4]>>[nH:1]1[cH:5][cH:6][nH:4][cH:7]1'
        return InchiParser.mol_from_inchi(product_inchi), smates

    @classmethod
    def test(cls):
        product, smates = cls.get_product_and_temp()
        product.draw()
        # temp.draw()
        for reactants in RxnTemplateApply.product_to_reactants(smates, product):
            for reactant in reactants:
                reactant.draw()
                print(SmilesWriter.mol_to_smiles(reactant))
                print(InchiWriter.mol_to_inchi(reactant))


if __name__ == "__main__":
    CrkitTester.test()
