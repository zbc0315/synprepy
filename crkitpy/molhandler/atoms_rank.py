# -*- coding: utf-8 -*-
# @Time     : 2020/5/28 14:48
# @Author   : zbc
# @Email    : zhangbc0315@outlook.com
# @FileName : atoms_rank.py
# @Software : PyCharm

from typing import List, Dict

import networkx as nx

from crkitpy.molgraph.molecule import Molecule, Atom, Bond


class AtomsRank:

    @classmethod
    def rank_atoms(cls, molecule: Molecule, for_template: bool = False):
        atoms_score: Dict[Atom, float] = nx.pagerank(molecule.mol_graph,
                                                      nstart=cls._nstart(molecule), weight="bond_weight")
        sorted_atoms = sorted(atoms_score.keys(), key=lambda x: atoms_score[x])
        for n, atom in enumerate(sorted_atoms):
            atom.set_id(n)

    @classmethod
    def _nstart(cls, molecule: Molecule) -> Dict[Atom, float]:
        nstart = {}
        for atom in molecule.get_atoms():
            nstart[atom] = float(f'{atom.get_atomic_num()}.'
                                 f'{5+atom.charge}'
                                 f'{atom.get_explicit_valence()}'
                                 f'{atom.get_implicit_valence()}'
                                 f'{atom.get_hs_num()}'
                                 f'{1 if atom.is_aromatic else 0}'.replace('-', ''))
        return nstart

    @classmethod
    def _nstart_for_temp(cls, template: Molecule) -> Dict[Atom, float]:
        nstart = cls._nstart(template)

        invisible_to_visible = {}
        for bond in template.get_used_linked_bonds():
            if not bond.get_atoms()[0].is_visible():
                cls._add_max_value_to_dict(invisible_to_visible, bond.get_atoms()[0], bond.get_atoms()[1], nstart)
            elif not bond.get_atoms()[1].is_visible():
                cls._add_max_value_to_dict(invisible_to_visible, bond.get_atoms()[1], bond.get_atoms()[0], nstart)

        for atom in template.get_atoms():
            if not atom.is_visible:
                if atom in invisible_to_visible.keys():
                    nstart[atom] += invisible_to_visible[atom]/100000
        return nstart

    @classmethod
    def _add_max_value_to_dict(cls, dic: Dict[Atom, Atom], key_atom: Atom, value_atom: Atom, nstart: Dict[Atom, float]):
        if key_atom in dic.keys():
            if nstart[value_atom] > nstart[dic[key_atom]]:
                dic[key_atom] = value_atom
        else:
            dic[key_atom] = value_atom
