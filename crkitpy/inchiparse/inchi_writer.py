# -*- coding: utf-8 -*-
# @Time     : 2020/6/12 14:27
# @Author   : zbc
# @Email    : zhangbc0315@outlook.com
# @FileName : inchi_writer.py
# @Software : PyCharm

from rdkit.Chem import AllChem

from crkitpy.molgraph.molecule import Molecule
from crkitpy.smilesparse.smiles_writer import SmilesWriter


class InchiWriter:

    @classmethod
    def _mol_to_inchi_by_rdkit(cls, molecule: Molecule):
        smiles = SmilesWriter.mol_to_smiles(molecule)
        rd_mol = AllChem.MolFromSmiles(smiles)
        return AllChem.MolToInchi(rd_mol)

    @classmethod
    def mol_to_inchi(cls, molecule: Molecule):
        return cls._mol_to_inchi_by_rdkit(molecule)

    @classmethod
    def smiles_to_inchi(cls, smiles: str) -> str:
        rd_mol = AllChem.MolFromSmiles(smiles)
        return AllChem.MolToInchi(rd_mol)
