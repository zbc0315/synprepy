# -*- coding: utf-8 -*-
# @Time     : 2020/4/21 11:26
# @Author   : zbc
# @Email    : zhangbc0315@outlook.com
# @FileName : rxn_center_detector.py
# @Software : PyCharm

from typing import List, Iterator

from crkitpy.molgraph.reaction import Reaction
from crkitpy.molgraph.molecule import Molecule
from crkitpy.molgraph.atom import Atom
from crkitpy.const import ATOM_PROP, RXN_ROLE


class RxnCenterDetector:

    @classmethod
    def detect_rxn_center(cls, reaction: Reaction):
        cls._detect_reacted_atoms_in_rxn(reaction)

    @classmethod
    def _detect_reacted_atoms_in_rxn(cls, reaction: Reaction) -> None:
        cls._set_atoms_env_in_rxn(reaction)
        for reactant in reaction.get_reactants():
            for r_atom in reactant.get_atoms():
                if r_atom.map_num is None or r_atom.map_num == 0:
                    continue
                for p_atom in cls._get_atoms_by_map_num(reaction.get_products(), r_atom.map_num):
                    if not cls._has_same_prop(r_atom, p_atom) or not cls._has_same_env(r_atom, p_atom):
                        r_atom.set_property(ATOM_PROP.K_IS_REACTED, True)
                        p_atom.set_property(ATOM_PROP.K_IS_REACTED, True)
                        r_atom.set_property(ATOM_PROP.K_REACT_LEVEL, 0)
                        p_atom.set_property(ATOM_PROP.K_REACT_LEVEL, 0)
                    else:
                        r_atom.set_property(ATOM_PROP.K_IS_REACTED, False)
                        p_atom.set_property(ATOM_PROP.K_IS_REACTED, False)

    @classmethod
    def _has_same_prop(cls, r_atom: Atom, p_atom: Atom):
        if r_atom.charge != p_atom.charge:
            return False
        if r_atom.get_valence() != p_atom.get_valence():
            return False
        return True

    @classmethod
    def _has_same_env(cls, r_atom: Atom, p_atom: Atom) -> bool:
        # if r_atom.get_hs_num() != p_atom.get_hs_num():
        #    return False
        if r_atom.get_property(ATOM_PROP.K_ENV) != p_atom.get_property(ATOM_PROP.K_ENV):
            return False
        for nbr_r_atom in r_atom.get_neighbors():
            if nbr_r_atom.map_num is None or nbr_r_atom.map_num == 0:
                return False
            nbr_p_atom = cls._get_atom_by_map_num_from_atoms(nbr_r_atom.map_num, p_atom.get_neighbors())
            r_bond = r_atom.get_bond(nbr_r_atom)
            p_bond = p_atom.get_bond(nbr_p_atom)
            if r_bond.get_bond_type() != p_bond.get_bond_type():
                return False
        return True

    @classmethod
    def _get_atom_by_map_num_from_atoms(cls, map_num: int, atoms: Iterator[Atom]) -> Atom:
        for atom in atoms:
            if atom.map_num == map_num:
                return atom

    @classmethod
    def _set_atoms_env_in_rxn(cls, reaction: Reaction) -> None:
        # TODO 应更改为直接从atom.get_env()获得
        for reactant in reaction.get_reactants():
            cls._detect_atoms_env(reactant)
        for product in reaction.get_products():
            cls._detect_atoms_env(product)

    @staticmethod
    def _detect_atoms_env(molecule: Molecule) -> None:
        for atom in molecule.get_atoms():
            neighbors = atom.get_neighbors()
            neighbor_map_num = [neighbor.map_num if neighbor.map_num is not None else 0 for neighbor in neighbors]
            neighbor_map_num.sort()
            atom.set_property(ATOM_PROP.K_ENV, neighbor_map_num)

    @classmethod
    def _get_atoms_by_map_num(cls, molecules: List[Molecule], map_num: int) -> Iterator[Atom]:
        for molecule in molecules:
            for atom in molecule.get_atoms():
                if atom.map_num == map_num:
                    yield atom

    @staticmethod
    def _get_molecules_by_role(reaction: Reaction, rxn_role: str) -> List[Molecule]:
        if rxn_role == RXN_ROLE.REACTANT:
            return reaction.get_reactants()
        elif rxn_role == RXN_ROLE.PRODUCT:
            return reaction.get_products()
        elif rxn_role == RXN_ROLE.SPECTATOR:
            return reaction.get_spectators()
        else:
            return []

