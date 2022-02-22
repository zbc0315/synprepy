#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Time    : 2022/1/5 12:40
# @Author  : zhangbc0315@outlook.com
# @File    : rxn_handler.py
# @Software: PyCharm

from typing import Iterator

from rdkit.Chem.rdchem import Mol, EditableMol
from rdkit.Chem import AllChem

from chem_utils.chem_handler import MolHandler


class RxnHandler:

    # region ===== smiles smarts =====

    @classmethod
    def smarts_to_rxn(cls, smarts: str, remove_stereo: bool):
        rxn = AllChem.ReactionFromSmarts(smarts)
        new_rxn = AllChem.ChemicalReaction()
        for reactant in rxn.GetReactants():
            if remove_stereo:
                AllChem.RemoveStereochemistry(reactant)
            new_rxn.AddReactantTemplate(AllChem.MolFromSmiles(AllChem.MolToSmiles(reactant)))
        for product in rxn.GetProducts():
            if remove_stereo:
                AllChem.RemoveStereochemistry(product)
            new_rxn.AddProductTemplate(AllChem.MolFromSmiles(AllChem.MolToSmiles(product)))
        return new_rxn

    # endregion

    # region ===== remove mapped mols =====

    @classmethod
    def _is_mol_mapped(cls, mol: Mol) -> bool:
        for atom in mol.GetAtoms():
            if atom.GetAtomMapNum() != 0:
                return True
        return False

    @classmethod
    def _get_mapped_mols(cls, mols: [Mol]) -> Iterator[Mol]:
        for mol in mols:
            if cls._is_mol_mapped(mol):
                yield mol

    @classmethod
    def remove_unmapped_mols_in_rxn(cls, rxn: AllChem.ChemicalReaction) -> AllChem.ChemicalReaction:
        new_rxn = AllChem.ChemicalReaction()
        for mapped_reactant in cls._get_mapped_mols(rxn.GetReactants()):
            new_rxn.AddReactantTemplate(mapped_reactant)
        for mapped_product in cls._get_mapped_mols(rxn.GetProducts()):
            new_rxn.AddProductTemplate(mapped_product)
        return new_rxn

    # endregion

    # region ===== remove products same with reactants

    @classmethod
    def remove_products_same_with_reactants(cls, rxn: AllChem.ChemicalReaction) -> AllChem.ChemicalReaction:
        reactants_smileses = [MolHandler.get_smiles_without_map_num(mol) for mol in rxn.GetReactants()]
        new_rxn = AllChem.ChemicalReaction()
        has_product = False
        for product in rxn.GetProducts():
            p_smiles = MolHandler.get_smiles_without_map_num(product)
            if p_smiles in reactants_smileses:
                continue
            new_rxn.AddProductTemplate(product)
            has_product = True
        if has_product:
            for reactant in rxn.GetReactants():
                new_rxn.AddReactantTemplate(reactant)
        return new_rxn

    # endregion


if __name__ == "__main__":
    pass
