# -*- coding: utf-8 -*-
# @Time     : 2020/6/11 14:20
# @Author   : zbc
# @Email    : zhangbc0315@outlook.com
# @FileName : aromatic_handler.py
# @Software : PyCharm

from typing import List, Iterator, Tuple

from crkitpy.molgraph.molecule import Molecule, AtomSet, Atom, Bond
from crkitpy.molhandler.aromatic_cycle_detector import AromaticCycleDetector
from crkitpy.molhandler.molecule_checker import MoleculeChecker
from crkitpy.const import BOND_TYPE
from crkitpy.molhandler.molecule_hs_num_corrector import MoleculeHsNumCorrector
from crkitpy.errors.bond_error import *


class AromaticHandler:

    @classmethod
    def remove_aromatic(cls, molecule: Molecule, need_correct_hs_num: bool = False):
        """ 移除分子中的芳香键，改为单键或双键或三键
        此方法之前必须先运行atom.correct_implicit_hs_num

        :param molecule:
        :param need_correct_hs_num:
        :return:
        """
        if need_correct_hs_num:
            MoleculeHsNumCorrector.correct_hs_num(molecule)
        AromaticCycleDetector.detect_aromatic_cycles(molecule)
        # 逐个处理每个芳香环
        for aromatic_cycle in cls._get_aromatic_cycles(molecule):
            passed_bonds = []
            while True:
                # 逐个得到芳环中每一个键，优先得到已经不是芳香键的键
                for bond in cls._get_bonds(aromatic_cycle):
                    if bond in passed_bonds:
                        continue
                    if bond.get_bond_type() == BOND_TYPE.AROMATIC:
                        bond.set_bond_type(BOND_TYPE.SINGLE)
                    passed_bonds.append(bond)
                    for atom, nbr_bond in cls._get_neighbor_cycle_bonds(bond, aromatic_cycle.get_bonds()):
                        if nbr_bond.get_bond_type() != BOND_TYPE.AROMATIC:
                            continue
                        cls._reset_bond_type(bond, nbr_bond, atom)
                if len(passed_bonds) == len(aromatic_cycle.get_bonds()):
                    break
        MoleculeChecker.check_aromatic_remove(molecule)

    @classmethod
    def _get_aromatic_cycles(cls, molecule: Molecule) -> Iterator[AtomSet]:
        # TODO 考虑到五元杂环中单双键的复杂情况，此处先处理偶数原子的芳环，是否能得到预期结果还有待证实
        odd_cycles = []
        for cycle in molecule.get_aromatic_cycles():
            if len(cycle.get_atoms()) % 2 == 0:
                yield cycle
            else:
                odd_cycles.append(cycle)
        for cycle in odd_cycles:
            yield cycle

    @classmethod
    def _reset_bond_type(cls, bond: Bond, nbr_bond: Bond, atom: Atom) -> None:
        """ 重置键的类型

        :param bond: 已被重置的键
        :param nbr_bond: 待重置的键
        :param atom: bond与nbr_bond之间的atom
        :return:
        """
        # atom.get_molecule().draw()
        another_atom = nbr_bond.get_another_atom(atom)
        rest_max_valence = cls._get_rest_max_valence(another_atom, nbr_bond)
        other_bonds_weight, aromatic_num = cls._get_other_bonds_weight(atom, [bond, nbr_bond])
        bond_weight = atom.get_valence() \
                      - atom.get_hs_num() \
                      - other_bonds_weight \
                      - BOND_TYPE.BOND_WEIGHTS[bond.get_bond_type()]
        if bond_weight > rest_max_valence:
            if aromatic_num > 0 and bond_weight == rest_max_valence + 1:
                bond_weight = rest_max_valence
            else:
                # pass
                raise BondValenceError('expect lighter bond')
        # bond_weight = bond_weight if bond_weight > 0 else 1
        nbr_bond.set_bond_type(BOND_TYPE.BOND_WEIGHT_STR_TO_SYMBOLS[str(bond_weight)])

    @classmethod
    def _get_rest_max_valence(cls, atom: Atom, bond: Bond) -> int:
        return int(atom.get_valence() - atom.get_hs_num() - cls._get_other_bonds_weight(atom, [bond])[0])

    @classmethod
    def _get_other_bonds_weight(cls, atom: Atom, bonds: List[Bond]) -> Tuple[float, int]:
        bonds_weight = 0
        aromatic_num = 0
        for bond in atom.get_bonds():
            if any([a.get_symbol() == 'H' for a in bond.get_atoms()]):
                continue
            if bond not in bonds:
                if bond.get_bond_type() == BOND_TYPE.AROMATIC:
                    bonds_weight += 1
                    aromatic_num += 1
                else:
                    bonds_weight += BOND_TYPE.BOND_WEIGHTS[bond.get_bond_type()]
        return bonds_weight, aromatic_num

    @classmethod
    def _get_neighbor_cycle_bonds(cls, bond: Bond, cycle_bonds: List[Bond]) -> Iterator[Tuple[Atom, Bond]]:
        for atom in bond.get_atoms():
            for nbr_bond in atom.get_bonds():
                if nbr_bond != bond and nbr_bond in cycle_bonds:
                    yield atom, nbr_bond

    @classmethod
    def _get_bonds(cls, cycle: AtomSet) -> Iterator[Bond]:
        for bond in cycle.get_bonds():
            if bond.get_bond_type() != BOND_TYPE.AROMATIC:
                yield bond
        yield cycle.get_bonds()[0]
