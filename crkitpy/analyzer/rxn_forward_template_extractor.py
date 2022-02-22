# -*- coding: utf-8 -*-
# @Time     : 2020/4/21 14:25
# @Author   : zbc
# @Email    : zhangbc0315@outlook.com
# @FileName : rxn_forward_template_extractor.py
# @Software : PyCharm

import networkx as nx

from crkitpy.rxnhandler.rxn_center_detector import RxnCenterDetector
from crkitpy.molgraph.reaction import Reaction
from crkitpy.molgraph.molecule import Molecule
from crkitpy.molgraph.atom import Atom
from crkitpy.const import ATOM_PROP


class RxnForwardTemplateExtractor:

    @classmethod
    def detect_rxn_template(cls, reaction: Reaction) -> Reaction:
        RxnCenterDetector.detect_rxn_center(reaction)
        cls._detect_need_save_atoms_in_rxn(reaction)
        return cls._extract_rxn_template(reaction)

    @classmethod
    def _extract_rxn_template(cls, reaction: Reaction) -> Reaction:
        rxn_template = reaction.copy()
        cls._replace_unsave_atoms_by_rgroup_in_rxn(rxn_template)
        return rxn_template

    @classmethod
    def _replace_unsave_atoms_by_rgroup_in_rxn(cls, reaction: Reaction) -> None:
        for reactant in reaction.get_reactants():
            cls._replace_unsave_atoms_by_rgroup_in_mol(reactant)
        for product in reaction.get_products():
            cls._replace_unsave_atoms_by_rgroup_in_mol(product)

    @classmethod
    def _replace_unsave_atoms_by_rgroup_in_mol(cls, molecule: Molecule) -> None:

        for saved_atom in molecule.get_saved_atom():
            for neighbor_atom in saved_atom.get_neighbors():
                if neighbor_atom.get_property(ATOM_PROP.K_NEED_SAVE):
                    continue
                bond = molecule.get_bond(saved_atom, neighbor_atom)
                new_bond = bond.copy()
                r_group = Atom.instance_by_symbol('*')
                r_group.set_property(ATOM_PROP.K_NEED_SAVE, True)
                molecule.add_atom(r_group)
                molecule.add_bond(new_bond, saved_atom, r_group)
                molecule.remove_bond(bond)

        cls._drop_unsaved_atoms_in_mol(molecule)

    @staticmethod
    def _drop_unsaved_atoms_in_mol(molecule: Molecule) -> None:
        atoms = molecule.get_atoms().copy()
        for atom in atoms:
            if not atom.get_property(ATOM_PROP.K_NEED_SAVE) or atom.get_property(ATOM_PROP.K_NEED_SAVE) is None:
                molecule.remove_atom(atom)

    @classmethod
    def _detect_need_save_atoms_in_rxn(cls, reaction: Reaction) -> None:
        cls._save_neighbor_atoms_in_rxn(reaction)
        cls._save_atoms_between_reacted_atoms_in_rxn(reaction)
        cls._save_atoms_in_product_map_to_reactant_saved_atoms(reaction)

    @classmethod
    def _save_neighbor_atoms_in_rxn(cls, reaction: Reaction) -> None:
        for reactant in reaction.get_reactants():
            cls._save_neighbor_atoms_in_mol(reactant)
#        for product in reaction.get_products():
#            cls._save_neighbor_atoms_in_mol(product)

    @classmethod
    def _save_atoms_between_reacted_atoms_in_rxn(cls, reaction: Reaction) -> None:
        for reactant in reaction.get_reactants():
            cls._save_atoms_between_reacted_atoms_in_mol(reactant)
#        for product in reaction.get_products():
#            cls._save_atoms_between_reacted_atoms_in_mol(product)

    @staticmethod
    def _save_atoms_in_product_map_to_reactant_saved_atoms(reaction: Reaction) -> None:
        for reactant in reaction.get_reactants():
            for r_atom in reactant.get_atoms():
                if r_atom.get_property(ATOM_PROP.K_NEED_SAVE):

                    for product in reaction.get_products():
                        for p_atom in product.get_atoms_by_map_num(r_atom.map_num):
                            p_atom.set_property(ATOM_PROP.K_NEED_SAVE, True)

    @staticmethod
    def _save_atoms_between_reacted_atoms_in_mol(molecule: Molecule) -> None:
        reacted_atoms = list(molecule.get_reacted_atom())
        for i in range(len(reacted_atoms)):
            for j in range(i+1, len(reacted_atoms)):
                atoms_in_path = nx.dijkstra_path(molecule.mol_graph, reacted_atoms[i], reacted_atoms[j])
                for atom in atoms_in_path:
                    atom.set_property(ATOM_PROP.K_NEED_SAVE, True)

    @staticmethod
    def _save_neighbor_atoms_in_mol(molecule: Molecule) -> None:
        for atom in molecule.get_atoms():
            if not atom.get_property(ATOM_PROP.K_IS_REACTED):
                continue
            atom.set_property(ATOM_PROP.K_NEED_SAVE, True)
            for neighbor in atom.get_neighbors():
                neighbor.set_property(ATOM_PROP.K_NEED_SAVE, True)
