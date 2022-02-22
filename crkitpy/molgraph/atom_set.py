# -*- coding: utf-8 -*-
# @Time     : 2020/5/29 16:18
# @Author   : zbc
# @Email    : zhangbc0315@outlook.com
# @FileName : atom_set.py
# @Software : PyCharm

from typing import List

from crkitpy.molgraph.atom import Atom
from crkitpy.molgraph.bond import Bond


class AtomSet:

    def __init__(self, atoms: List[Atom] = None):
        if atoms is None:
            self._atoms: List[Atom] = []
        else:
            self._atoms = atoms
        self._bonds = []

    def get_atoms(self) -> List[Atom]:
        return self._atoms

    def add_atom(self, atom: Atom) -> None:
        if atom not in self._atoms:
            self._atoms.append(atom)

    def add_bond(self, bond: Bond) -> None:
        if bond not in self._bonds:
            self._bonds.append(bond)

    def get_bonds(self) -> List[Bond]:
        return self._bonds
