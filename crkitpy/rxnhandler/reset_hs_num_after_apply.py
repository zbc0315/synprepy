#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Time    : 2021/9/17 10:03
# @Author  : zhangbc0315@outlook.com
# @File    : reset_hs_num_after_apply.py
# @Software: PyCharm

from crkitpy.const import ATOM_PROP, BOND_TYPE
from crkitpy.const.bond_type import BOND_WEIGHTS_TO_TYPE
from crkitpy.molhandler.cycle_detector import CycleDetector
from crkitpy.errors.bond_error import BondError


class ResetHsNumAfterApply:

    # ================
    # -- 获得键的贡献 --

    @classmethod
    def get_bonds_contributes(cls, atom, except_bonds=None):
        res = 0
        except_bonds = [] if except_bonds is None else except_bonds
        aromatic_bonds = []
        for bond in atom.get_bonds():
            if bond in except_bonds:
                continue
            if bond.get_bond_type() == BOND_TYPE.AROMATIC:
                aromatic_bonds.append(bond)
                continue
            res += bond.get_bond_contribute()
        if len(aromatic_bonds) == 1:
            raise BondError("应用完模板，修正H原子数时发现，某反应中心原子只连接了1个芳香键")
        elif len(aromatic_bonds) == 2:
            if atom.get_atomic_num() == 7:
                res += 2
            else:
                res += 3
        elif len(aromatic_bonds) == 3:
            res += 3
        elif len(aromatic_bonds) > 3:
            raise BondError(f"应用完模板，修正H原子数时发现，某反应中心原子连接了{len(aromatic_bonds)}个芳香键")
        return res

    # =====================
    # -- 查询环中的原子或键 --

    @classmethod
    def _get_all_bonds_in_cycle_by_atoms(cls, atoms):
        bonds = set()
        for atom in atoms:
            for bond in atom.get_bonds():
                if bond.get_another_atom(atom) in atoms:
                    bonds.add(bond)
        return list(bonds)

    # =============================
    # -- 将不在芳环内的原子消除芳香性 --

    @classmethod
    def _set_atoms_not_in_aro_cycle_no_aro(cls, mol, atoms_in_cycles, bonds_in_cycles):
        all_atoms_in_cycles = set()
        for atoms in all_atoms_in_cycles:
            all_atoms_in_cycles = all_atoms_in_cycles.union(atoms)
        for atom in mol.get_atoms():
            if not atom.is_aromatic:
                continue
            if atom not in all_atoms_in_cycles:
                atom.is_aromatic = False

    # =======================
    # -- 消除一些环内的芳香性 --

    @classmethod
    def _set_some_cycle_no_aromatic(cls, bonds_in_cycles):
        for bonds in bonds_in_cycles:
            if cls._bonds_should_set_no_aromatic(bonds, bonds_in_cycles):
                cls._set_cycle_no_aromatic(bonds, bonds_in_cycles)

    @classmethod
    def _set_cycle_no_aromatic(cls, bonds, bonds_in_cycles):
        certain_bonds = cls._get_certain_bonds(bonds, bonds_in_cycles)
        cls._set_some_easy_bond_in_cycle_no_aromatic(bonds, certain_bonds)
        num_need_change = len(bonds) - len(certain_bonds)
        num = 0
        while num_need_change > 0:
            num += 1
            if num >= 100:
                raise ValueError("尝试将非芳环中的芳香键转化为常规键，查找太多次")
            certain_atoms = cls._get_certain_atoms(certain_bonds)
            bond = cls._get_bond_nbr_certain(bonds, certain_bonds, certain_atoms)
            if bond is None:
                break
            for atom in bond.get_atoms():
                if atom not in certain_atoms or atom.get_property(ATOM_PROP.K_REACT_CENTER):
                    continue
                other_bond_level = cls.get_bonds_contributes(atom, [bond])
                bond_level = atom.get_valence() - other_bond_level - atom.charge - atom.get_hs_num()
                bond.set_bond_type(BOND_WEIGHTS_TO_TYPE[bond_level])
                certain_bonds.append(bond)
                num_need_change -= 1
                break

    @classmethod
    def _set_some_easy_bond_in_cycle_no_aromatic(cls, bonds, certain_bonds):
        atoms = set()
        for bond in bonds:
            atoms = atoms.union(bond.get_atoms())
        num = 0
        while True:
            num += 1
            if num >= 10:
                raise ValueError("尝试将非芳环中芳香键转化为普通键，迭代次数超过10次")
            had_changed = False
            for atom in atoms:
                if atom.get_property(ATOM_PROP.K_REACT_CENTER) is None:
                    atom_certain_bonds = []
                    atom_uncertain_bonds = []
                    for bond in atom.get_bonds():
                        if bond not in bonds or bond in certain_bonds:
                            atom_certain_bonds.append(bond)
                        else:
                            atom_uncertain_bonds.append(bond)
                    certain_bonds_contribute = cls.get_bonds_contributes(atom, atom_uncertain_bonds)
                    rest_level = atom.get_valence() - certain_bonds_contribute - atom.get_hs_num() - atom.charge
                    if len(atom_uncertain_bonds) >= rest_level > 0:
                        had_changed = True
                        for bond in atom_uncertain_bonds:
                            bond.set_bond_type(BOND_TYPE.SINGLE)
                            certain_bonds.append(bond)
            if not had_changed:
                break

    @classmethod
    def _get_certain_atoms(cls, certain_bonds):
        certain_atoms = set()
        for c_bond in certain_bonds:
            certain_atoms = certain_atoms.union(set(c_bond.get_atoms()))
        return certain_atoms

    @classmethod
    def _get_bond_nbr_certain(cls, bonds, certain_bonds, certain_atoms):
        for bond in bonds:
            if bond in certain_bonds:
                continue
            if any([atom in certain_atoms and atom.get_property(ATOM_PROP.K_REACT_CENTER) is None for atom in bond.get_atoms()]):
                return bond

    @classmethod
    def _get_certain_bonds(cls, bonds, bonds_in_cycle):
        res = []
        for bond in bonds:
            if bond.get_bond_type() != BOND_TYPE.AROMATIC:
                res.append(bond)
            elif cls._bond_in_aromatic_cycle(bond, bonds_in_cycle):
                res.append(bond)
        return res

    @classmethod
    def _bonds_should_set_no_aromatic(cls, bonds, bonds_in_cycles):
        if any([bond.get_bond_type() != BOND_TYPE.AROMATIC for bond in bonds]):
            for bond in bonds:
                if bond.get_bond_type() == BOND_TYPE.AROMATIC and not cls._bond_in_aromatic_cycle(bond, bonds_in_cycles):
                    return True
        return False

    @classmethod
    def _bond_in_aromatic_cycle(cls, bond, bonds_in_cycles):
        for bonds in bonds_in_cycles:
            if bond in bonds and all([bond.get_bond_type() == BOND_TYPE.AROMATIC for bond in bonds]):
                return True
        return False

    # ==============================
    # -- 将不在环中芳香键设置为非芳香键 --

    @classmethod
    def _set_bonds_not_in_cycle_no_aromatic(cls, mol, bonds_in_cycles):
        all_bonds_in_cycles = set()
        for bonds in bonds_in_cycles:
            all_bonds_in_cycles = all_bonds_in_cycles.union(set(bonds))
        for bond in mol.get_bonds():
            if bond in all_bonds_in_cycles:
                continue
            if bond.get_bond_type() != BOND_TYPE.AROMATIC:
                continue
            bond_level = None
            for atom in bond.get_atoms():
                if atom.get_property(ATOM_PROP.K_REACT_CENTER):
                    continue
                other_bond_level = cls.get_bonds_contributes(atom, [bond])
                bond_level = atom.get_valence() - other_bond_level - atom.charge - atom.get_hs_num()
                break
            bond.set_bond_type(BOND_WEIGHTS_TO_TYPE[bond_level])

    # ==============================
    # -- 将芳香环中的所有键设置为芳香键 --

    @classmethod
    def _cycle_is_aromatic(cls, atoms):
        if all([atom.is_aromatic for atom in atoms]):
            return True

    @classmethod
    def _set_bonds_aromatic(cls, bonds):
        for bond in bonds:
            bond.set_bond_type(BOND_TYPE.AROMATIC)

    @classmethod
    def _reset_aromatic_cycles(cls, mol):
        all_bonds = []
        all_atoms = []
        for cycle in CycleDetector.get_cycles(mol):
            bonds = cls._get_all_bonds_in_cycle_by_atoms(cycle)
            all_bonds.append(bonds)
            all_atoms.append(cycle)
            if cls._cycle_is_aromatic(cycle):
                cls._set_bonds_aromatic(bonds)
        return all_bonds, all_atoms

    # =======================

    @classmethod
    def reset_hs_num(cls, reactants):
        """ 重置H原子数
        1. 如果环中所有原子都是芳香，将所有键设置为芳香键
        2. 如果有芳香键不在环中，将它设置为非芳香键
        3. 如果有些环，只有部分独占键是芳香键，消除环的芳香性

        """
        for reactant in reactants:
            # reactant.draw()
            bonds_in_cycles, atoms_in_cycles = cls._reset_aromatic_cycles(reactant)
            # reactant.draw()
            cls._set_bonds_not_in_cycle_no_aromatic(reactant, bonds_in_cycles)
            # reactant.draw()
            cls._set_some_cycle_no_aromatic(bonds_in_cycles)
            # reactant.draw()
            cls._set_atoms_not_in_aro_cycle_no_aro(reactant, atoms_in_cycles, bonds_in_cycles)
            # reactant.draw()
            for atom in reactant.get_atoms():
                if atom.get_property(ATOM_PROP.K_REACT_CENTER) is None:
                    continue
                bonds_contribute = cls.get_bonds_contributes(atom)
                hs_num = atom.get_valence() - bonds_contribute - atom.charge
                atom.implicit_hs_num = hs_num - atom.explicit_hs_num


if __name__ == "__main__":
    pass
