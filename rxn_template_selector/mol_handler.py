# -*- coding: utf-8 -*-
# @Time     : 2021/4/25 11:05
# @Author   : zbc
# @Email    : zhangbc0315@outlook.com
# @FileName : mol_handler.py
# @Software : PyCharm

from rdkit.Chem import AllChem


class MolHandler:

    @classmethod
    def frag_by_aids(cls, mol, aids: [int]):
        edit_mol = AllChem.EditableMol(mol)
        for atom in mol.GetAtoms():
            if atom.GetIdx() not in aids:
                edit_mol.RemoveAtom(atom.GetIdx())
        mol = edit_mol.GetMol()
        AllChem.SanitizeMol(mol)
        return mol

    @classmethod
    def smarts_of_frag_by_aids(cls, mol: AllChem.Mol, aids: [int]):
        # return AllChem.MolToSmiles(cls.frag_by_aids(mol, aids))
        smi = AllChem.MolFragmentToSmarts(mol, aids)
        # mol = AllChem.MolFromSmiles(smi)
        # AllChem.SanitizeMol(mol)
        # smi = AllChem.MolToSmiles(mol)
        return smi


if __name__ == '__main__':
    m = AllChem.MolFromSmiles("[CH3:8][C:1]1=[CH:6][CH:5]=[CH:4][CH:3]=[CH:2]1")
    print(MolHandler.smarts_of_frag_by_aids(m, [0, 1]))
