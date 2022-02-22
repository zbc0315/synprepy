# -*- coding: utf-8 -*-
# @Time     : 2020/6/10 12:47
# @Author   : zbc
# @Email    : zhangbc0315@outlook.com
# @FileName : molecule_checker.py
# @Software : PyCharm

from crkitpy.molgraph.molecule import Molecule, Atom
from crkitpy.const import BOND_TYPE, ATOM_PROP
from crkitpy.errors.atom_error import *
from crkitpy.errors.bond_error import *


class MoleculeChecker:

    # ===========================
    # -- check aromatic remove --

    @classmethod
    def check_aromatic_remove(cls, molecule: Molecule):
        for bond in molecule.get_bonds():
            if bond.get_bond_type() == BOND_TYPE.AROMATIC:
                bond.set_bond_type(BOND_TYPE.SINGLE)  # TODO 是否可以直接改成单键而不报错呢？？？
                # raise BondTypeError('unexpect aromatic bond')

    # ============================
    # -- check molecule valence --

    @classmethod
    def check_valence(cls, atom: Atom):
        if atom.get_valence() < atom.get_bonds_weight_without_h():
            # TODO 此处不应该直接取最大Valence
            if atom.get_property(ATOM_PROP.K_MAX_VALENCE) > atom.get_bonds_weight_without_h():
                atom.set_implicit_valence(atom.get_property(ATOM_PROP.K_MAX_VALENCE) - atom.get_explicit_valence())
            else:
                raise AtomValenceError('too much bond')

    @classmethod
    def check_molecule(cls, molecule: Molecule, need_correct_hs_num: True):
        for atom in molecule.get_atoms():
            if need_correct_hs_num:
                atom.correct_implicit_hs_num()
            cls.check_valence(atom)
