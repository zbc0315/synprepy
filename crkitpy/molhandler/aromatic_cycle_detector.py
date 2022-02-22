# -*- coding: utf-8 -*-
# @Time     : 2020/5/29 15:56
# @Author   : zbc
# @Email    : zhangbc0315@outlook.com
# @FileName : aromatic_cycle_detector.py
# @Software : PyCharm

from crkitpy.molhandler.cycle_detector import CycleDetector
from crkitpy.molgraph.molecule import Molecule, Atom, AtomSet


class AromaticCycleDetector:

    @classmethod
    def detect_aromatic_cycles(cls, molecule: Molecule):
        cycles = CycleDetector.get_cycles(molecule, cls._is_aromatic_atom)
        for cycle in cycles:
            cycle = AtomSet(cycle)
            molecule.add_aromatic_cycle(cycle)

    @classmethod
    def _is_aromatic_atom(cls, atom):
        return atom.is_aromatic
