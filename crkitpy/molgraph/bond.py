# -*- coding: utf-8 -*-
# @Time     : 2020/4/20 16:22
# @Author   : zbc
# @Email    : zhangbc0315@outlook.com
# @FileName : bond.py
# @Software : PyCharm

from typing import List, Tuple

from enum import Enum

from crkitpy.molgraph.mol_component import MolComponent
from crkitpy.const import BOND_TYPE
from crkitpy.const.bond_type import BOND_WEIGHTS


class BondType(Enum):
    UNSPECIFIC = 0
    SINGLE = 1
    DOUBLE = 2
    TRIPLE = 3
    AROMATIC = 1.5


class Bond(MolComponent):

    def __init__(self):
        super().__init__()
        self._atoms = []
        self.symbol = None
        self._bond_type = None
        self.bond_to_r = False

    def __str__(self):
        return str(self._bond_type)

    # ===============
    # -- bond type --

    def set_bond_type(self, bt: str) -> None:
        self._bond_type = bt
        self.symbol = BOND_TYPE.BOND_TO_SYMBOLS[self._bond_type]

    def get_bond_type(self) -> str:
        return self._bond_type

    def link_atoms(self, begin_atom, end_atom) -> None:
        self._molecule.add_bond(self, begin_atom, end_atom)

    def get_atoms(self) -> List:
        # TODO 检查个数，应该为两个
        return self._atoms

    def get_another_atom(self, atom):
        if atom not in self.get_atoms():
            raise ValueError('Cannot get another atom for atom not belong this bond')
        for _atom in self.get_atoms():
            if _atom != atom:
                return _atom

    def add_atom(self, atom) -> None:
        self._atoms.append(atom)

    def equal_in_type_and_atoms(self, other) -> bool:
        if self._bond_type != other.get_bond_type():
            return False
        return set(self.get_atoms()) == set(other.get_atoms())

    def is_link_atom_by_symbol(self, symbol: str) -> bool:
        for atom in self.get_atoms():
            if atom.get_symbol() == symbol:
                return True
        return False

    def get_bond_contribute(self):
        return BOND_WEIGHTS[self._bond_type]

    def copy(self, save_parent=False):
        bond = Bond()
        bond.set_bond_type(self._bond_type)
        bond.symbol = self.symbol
        if save_parent:
            bond.set_parent(self)
        return bond

    @staticmethod
    def instance_by_bond_type(bond_type: str):
        bond = Bond()
        bond.set_bond_type(bond_type)
        bond.symbol = BOND_TYPE.BOND_TO_SYMBOLS[bond_type]
        return bond
