# -*- coding: utf-8 -*-
# @Time     : 2020/6/24 15:03
# @Author   : zbc
# @Email    : zhangbc0315@outlook.com
# @FileName : inchi_parser.py
# @Software : PyCharm

from rdkit.Chem import AllChem

from crkitpy.smilesparse.rdkit_smiles_parser import RdkitSmilesParser


class InchiParser:

    @classmethod
    def mol_from_inchi(cls, inchi: str):
        rdmol = AllChem.MolFromInchi(inchi)
        return RdkitSmilesParser.rdmol_to_mol(rdmol)
