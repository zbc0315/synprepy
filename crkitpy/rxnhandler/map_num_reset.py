# -*- coding: utf-8 -*-
# @Time     : 2020/5/28 16:31
# @Author   : zbc
# @Email    : zhangbc0315@outlook.com
# @FileName : map_num_reset.py
# @Software : PyCharm

from typing import List

from crkitpy.molgraph.reaction import Reaction, Molecule
from crkitpy.molhandler.atoms_rank import AtomsRank
from crkitpy.errors.rxn_error import MoleculeNumError


class MapNumReset:

    # ============================
    # -- reset reaction map num --

    @classmethod
    def reset_rxn_map_num(cls, reaction: Reaction):
        pass

    @classmethod
    def _sort_map_num_by_id(cls, reaction: Reaction) -> List[int]:
        map_num_to_id = {}
        for product in reaction.get_products():
            for atom in product.get_atoms():
                if atom.map_num is not None and atom.map_num > 0:
                    if atom.map_num not in map_num_to_id.keys():
                        map_num_to_id[atom.map_num] = atom.get_id()
                    elif atom.get_id() < map_num_to_id[atom.map_num]:
                        map_num_to_id[atom.map_num] = atom.get_id()
        return sorted(map_num_to_id.keys(), key=lambda x: map_num_to_id[x])

    # =====================================
    # -- reset reaction template map num --

    @classmethod
    def reset_rxn_template_map_num(cls, rxn_template: Reaction, need_rank: bool = True):
        cls._check_is_rxn_template(rxn_template)
        if need_rank:
            AtomsRank.rank_atoms(rxn_template.get_reactants()[0], for_template=True)
        for r_atom in rxn_template.get_reactants()[0].get_atoms():
            if r_atom.map_num is not None and r_atom.map_num > 0:
                p_atom = rxn_template.get_products()[0].get_one_atom_by_map_num(map_num=r_atom.map_num)
                p_atom.set_id(r_atom.get_id())
            r_atom.map_num = r_atom.get_id() + 1
        for p_atom in rxn_template.get_products()[0].get_atoms():
            p_atom.map_num = p_atom.get_id() + 1

    @classmethod
    def _check_is_rxn_template(cls, rxn_template: Reaction):
        r_num = len(rxn_template.get_reactants())
        p_num = len(rxn_template.get_products())
        if r_num != 1 or p_num != 1:
            raise MoleculeNumError(f'except one reactant and one product in reaction template, '
                                   f'but get {r_num} reactant and {p_num} product')

    # ==================================
    # -- reset molecule map num to id --

    @classmethod
    def reset_molecule_map_num_to_id(cls, molecule: Molecule, need_rank: bool = True) -> None:
        if need_rank:
            AtomsRank.rank_atoms(molecule, for_template=True)
        for atom in molecule.get_atoms():
            atom.map_num = atom.get_id() + 1

    # ============================
    # -- reset molecule map num --

    @classmethod
    def reset_molecule_map_num(cls, molecule: Molecule, need_rank: bool = True) -> None:
        if need_rank:
            AtomsRank.rank_atoms(molecule)
        sorted_map_num = cls._sort_mol_map_num_by_id(molecule)
        atom_to_map_num = {}
        for n, map_num in enumerate(sorted_map_num):
            for atom in molecule.get_atoms_by_map_num(map_num):
                atom_to_map_num[atom] = n + 1
        for atom, map_num in atom_to_map_num.items():
            atom.map_num = map_num

    @classmethod
    def _sort_mol_map_num_by_id(cls, molecule: Molecule) -> List[int]:
        map_num_to_id = {}
        for atom in molecule.get_atoms():
            if atom.map_num is not None and atom.map_num > 0:
                if atom.map_num not in map_num_to_id.keys() or atom.get_id() < map_num_to_id[atom.map_num]:
                    map_num_to_id[atom.map_num] = atom.get_id()
        return sorted(map_num_to_id.keys(), key=lambda x: map_num_to_id[x])

