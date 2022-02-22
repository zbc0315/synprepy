# -*- coding: utf-8 -*-
# @Time     : 2021/6/17 14:57
# @Author   : zbc
# @Email    : zhangbc0315@outlook.com
# @FileName : rxn_mapper.py
# @Software : PyCharm

from rdkit.Chem import AllChem


class RxnMapper:

    @classmethod
    def _is_mol_mapped(cls, mol) -> bool:
        for atom in mol.GetAtoms():
            if atom.GetAtomMapNum() is not None and atom.GetAtomMapNum() != 0:
                return True
        return False

    @classmethod
    def _remove_map_num(cls, mol):
        for atom in mol.GetAtoms():
            atom.SetAtomMapNum(0)

    @classmethod
    def _clean_unmapped_mols(cls, mols) -> ([str], [str]):
        mapped = []
        unmapped = []
        for mol in mols:
            if cls._is_mol_mapped(mol):
                cls._remove_map_num(mol)
                mapped.append(mol)
            else:
                unmapped.append(mol)
        return mapped, unmapped

    @classmethod
    def get_mapped_reactants(cls, rxn_smi: str):
        rxn = AllChem.ReactionFromSmarts(rxn_smi)
        mapped_rs, unmapped_rs = cls._clean_unmapped_mols(rxn.GetReactants())
        return [AllChem.MolFromSmiles(AllChem.MolToSmiles(m)) for m in mapped_rs], [AllChem.MolToInchi(p) for p in rxn.GetProducts()]


if __name__ == '__main__':
    pass
