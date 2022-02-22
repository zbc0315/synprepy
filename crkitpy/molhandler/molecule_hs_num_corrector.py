# -*- coding: utf-8 -*-
# @Time     : 2020/6/11 15:53
# @Author   : zbc
# @Email    : zhangbc0315@outlook.com
# @FileName : molecule_hs_num_corrector.py
# @Software : PyCharm

from crkitpy.molgraph.molecule import Molecule, Atom


class MoleculeHsNumCorrector:

    @classmethod
    def correct_hs_num(cls, molecule: Molecule) -> None:
        for atom in molecule.get_atoms():
            atom.correct_implicit_hs_num()
