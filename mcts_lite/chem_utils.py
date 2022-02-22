#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Time    : 2021/9/13 11:39
# @Author  : zhangbc0315@outlook.com
# @File    : chem_utils.py
# @Software: PyCharm

from rdkit.Chem import AllChem


class ChemUtils:

    @classmethod
    def remove_inchi_stereo(cls, inchi):
        mol = AllChem.MolFromInchi(inchi)
        AllChem.RemoveStereochemistry(mol)
        return AllChem.MolToInchi(mol)


if __name__ == "__main__":
    pass
