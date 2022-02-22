# -*- coding: utf-8 -*-
# @Time     : 2020/7/14 16:18
# @Author   : zbc
# @Email    : zhangbc0315@outlook.com
# @FileName : calc_atoms_coordinate.py
# @Software : PyCharm

from typing import Iterator, List

from crkitpy.molgraph.molecule import Molecule
from crkitpy.molgraph.atom import Atom, Hybrid
from crkitpy.smilesparse.smiles_parser import SmilesParser
from crkitpy.molhandler.cycle_detector import CycleDetector
from crkitpy.const import ATOM_PROP


class CalcAtomsCoordinate:

    def __init__(self):
        passed_atoms = {}

    @classmethod
    def calc_2d_coordinate(cls, molecule: Molecule):
        cycles = CycleDetector.get_cycles(molecule)
        first_atom = molecule.get_atoms()[0]
        cls._traverse_atoms_and_set_alter_angles(first_atom, cycles)
        print(1)

    @classmethod
    def _get_owning_cycle(cls, atom: Atom, cycles: List[List[Atom]]):
        pass

    @classmethod
    def _traverse_atoms_and_set_alter_angles(cls, first_atom: Atom, cycles: List[List[Atom]]):
        """ 从第一个原子开始，遍历分子中每一个原子，设置原子间的前后关系, 设置备选键角

        :param first_atom: 第一个原子
        :return:
        """
        current_atoms = [first_atom]
        while True:
            next_atoms = []
            for current_atom in current_atoms:
                for next_atom in current_atom.get_neighbors():
                    if next_atom.get_property(ATOM_PROP.K_NEXT) is not None \
                            or next_atom.get_property(ATOM_PROP.K_PREVIOUS) is not None:
                        continue
                    next_atom.set_property(ATOM_PROP.K_PREVIOUS, current_atom)
                    if current_atom.get_property(ATOM_PROP.K_NEXT) is None:
                        current_atom.set_property(ATOM_PROP.K_NEXT, [next_atom])
                    else:
                        current_atom.get_property(ATOM_PROP.K_NEXT).append(next_atom)
                    next_atoms.append(next_atom)
            if len(next_atoms) == 0:
                break
            current_atoms = next_atoms.copy()


if __name__ == '__main__':
    smi = '[CH3:25][CH:21]([NH:20][P:19](=[O:26])([O:18][CH2:17][CH:14]1[O:13][C:8]([C:9]#[N:10])([CH:11]([OH:12])[CH:15]1[OH:16])[c:7]1[cH:34][cH:35][c:36]2[c:2]([NH2:1])[nH:3][cH:4][nH:5][n:6]12)[O:27][c:28]1[cH:29][cH:30][cH:31][cH:32][cH:33]1)[C:22]([OH:24])=[O:23]'
    mol = SmilesParser.parse_to_mol_by_lex(smi)
    CalcAtomsCoordinate.calc_2d_coordinate(mol)
