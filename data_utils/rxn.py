#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Time    : 2022/2/11 14:58
# @Author  : zhangbc0315@outlook.com
# @File    : rxn.py
# @Software: PyCharm

from typing import Optional

from data_utils.mol_data import Mol, MolData


class Rxn:

    def __init__(self, rxn_code: str,
                 rxn_mols_code: str,
                 reactants_codes: [str],
                 product_code: str,
                 catalysts_codes: [str],
                 solvents_codes: [str],
                 yield_product: float,
                 rxn_smiles: str,
                 create_year: int):
        self._rxn_code: str = rxn_code
        self._rxn_mols_code: str = rxn_mols_code
        self._reactants_codes: [str] = reactants_codes
        self._product_code: str = product_code
        self._catalysts_codes: [str] = catalysts_codes
        self._solvents_codes: [str] = solvents_codes
        self._yield_product: float = yield_product
        self._rxn_smiles: str = rxn_smiles
        self._create_year: int = create_year

        self._reactants: [Mol] = []
        self._product: Mol = Optional[None]
        self._catalysts: [Mol] = []
        self._solvents: [Mol] = []

        self._rct_code: str = Optional[None]
        self._ret_code: str = Optional[None]

    def init_rxn(self, mol_data: MolData):
        for mol_code in self._reactants_codes:
            self._reactants.append(mol_data.get_mol_by_mol_code(mol_code))
        self._product = mol_data.get_mol_by_mol_code(self._product_code)
        for mol_code in self._catalysts_codes:
            self._catalysts.append(mol_data.get_mol_by_mol_code(mol_code))
        for mol_code in self._solvents_codes:
            self._solvents.append(mol_data.get_mol_by_mol_code(mol_code))

    # region ===== rxn_code =====
    @property
    def rxn_code(self):
        return self._rxn_code

    @rxn_code.setter
    def rxn_code(self, value):
        self._rxn_code = value
    # endregion

    # region ===== rxn_smiles =====
    @property
    def rxn_smiles(self):
        return self._rxn_smiles

    @rxn_smiles.setter
    def rxn_smiles(self, value):
        self._rxn_smiles = value
    # endregion

    # region ===== yield_product =====
    @property
    def yield_product(self):
        return self._yield_product

    @yield_product.setter
    def yield_product(self, value):
        self._yield_product = value
    # endregion

    # region ===== create_year =====
    @property
    def create_year(self):
        return self._create_year

    @create_year.setter
    def create_year(self, value):
        self._create_year = value
    # endregion

    # region ===== rxn_mols_code =====
    @property
    def rxn_mols_code(self):
        return self._rxn_mols_code

    @rxn_mols_code.setter
    def rxn_mols_code(self, value):
        self._rxn_mols_code = value
    # endregion

    # region ===== mol =====
    @property
    def reactants(self) -> [Mol]:
        return self._reactants

    @property
    def product(self) -> Mol:
        return self._product

    @property
    def catalysts(self) -> [Mol]:
        return self._catalysts

    @property
    def solvents(self) -> [Mol]:
        return self._solvents

    @property
    def catalysts_codes(self) -> [str]:
        return self._catalysts_codes

    @property
    def solvents_codes(self) -> [str]:
        return self._solvents_codes
    # endregion

    # region ===== rct_code =====
    @property
    def rct_code(self):
        return self._rct_code

    @rct_code.setter
    def rct_code(self, value):
        self._rct_code = value
    # endregion

    # region ===== ret_code =====
    @property
    def ret_code(self):
        return self._ret_code

    @ret_code.setter
    def ret_code(self, value):
        self._ret_code = value
    # endregion


if __name__ == "__main__":
    pass
