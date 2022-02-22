#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Time    : 2022/1/5 11:21
# @Author  : zhangbc0315@outlook.com
# @File    : mol_handler.py
# @Software: PyCharm

from typing import Iterator, Tuple, Dict

from rdkit.Chem.rdchem import Mol, EditableMol, Atom
from rdkit.Chem import AllChem


class MolHandler:

    # region ===== get mol fingerprint =====

    @classmethod
    def get_mol_fp_by_rdmol(cls, rdmol: Mol):
        return AllChem.GetMorganFingerprintAsBitVect(rdmol)

    @classmethod
    def get_mol_fp_by_smiles(cls, smiles: str):
        return cls.get_mol_fp_by_rdmol(AllChem.MolFromSmiles(smiles))

    # endregion

    # region ===== get smiles =====

    @classmethod
    def get_smiles_without_map_num(cls, mol: Mol) -> str:
        copy_mol = AllChem.MolFromSmiles(AllChem.MolToSmiles(mol))
        for atom in copy_mol.GetAtoms():
            atom.SetAtomMapNum(0)
        return AllChem.MolToSmiles(copy_mol)

    # endregion

    # region ===== canonicalize map num =====

    @classmethod
    def canonicalize_map_num(cls, mol: Mol) -> Dict[int, int]:
        map_num_old_to_new = {}
        new_map_num = 1
        ranked_atom_idxes = list(AllChem.CanonicalRankAtoms(mol))
        for idx in ranked_atom_idxes:
            atom = mol.GetAtomWithIdx(idx)
            old_map_num = atom.GetAtomMapNum()
            if old_map_num == 0:
                continue
            if old_map_num in map_num_old_to_new.keys():
                atom.SetAtomMapNum(map_num_old_to_new[old_map_num])
            else:
                atom.SetAtomMapNum(new_map_num)
                map_num_old_to_new[old_map_num] = new_map_num
                new_map_num += 1
        AllChem.CanonicalizeMol(mol)
        return map_num_old_to_new

    # endregion

    # region ===== edit mol =====

    @classmethod
    def merge_mols(cls, mols: [Mol]) -> Mol:
        res = None
        for n, mol in enumerate(mols):
            if n == 0:
                res = mol
            else:
                res = AllChem.CombineMols(res, mol)
        return res

    @classmethod
    def frag_idxes_to_mol(cls, mol: Mol, idxes: [int]):
        em = EditableMol(Mol())
        idxes_old_to_new = {}
        for i, idx in enumerate(idxes):
            em.AddAtom(mol.GetAtomWithIdx(idx))
            idxes_old_to_new[idx] = i

        for i, idx in enumerate(idxes):
            atom = mol.GetAtomWithIdx(idx)
            for bond in atom.GetBonds():
                if bond.GetBeginAtomIdx() == idx:
                    other_idx = bond.GetEndAtomIdx()
                else:
                    other_idx = bond.GetBeginAtomIdx()
                if other_idx < idx:
                    continue
                if other_idx in idxes_old_to_new.keys():
                    em.AddBond(idxes_old_to_new[idx], idxes_old_to_new[other_idx], bond.GetBondType())
        res = em.GetMol()
        res.ClearComputedProps()
        AllChem.GetSymmSSSR(res)
        res.UpdatePropertyCache(False)
        res._idxMap = idxes_old_to_new
        return res

    # endregion

    # region ===== query atom with map num =====

    @classmethod
    def contain_atom_with_map_num_in_mols(cls, mols: [Mol], map_num: int) -> bool:
        for _ in cls.query_atoms_with_map_num_in_mols(mols, map_num):
            return True
        return False

    @classmethod
    def query_atoms_with_map_num_in_mols(cls, mols: [Mol], map_num: int) -> Iterator[Atom]:
        for mol in mols:
            for atom in cls.query_atoms_with_map_num(mol, map_num):
                yield atom

    @classmethod
    def query_atoms_with_map_num(cls, mol: Mol, map_num: int) -> Iterator[Atom]:
        for atom in mol.GetAtoms():
            if atom.GetAtomMapNum() == map_num:
                yield atom

    # endregion


if __name__ == "__main__":
    pass
