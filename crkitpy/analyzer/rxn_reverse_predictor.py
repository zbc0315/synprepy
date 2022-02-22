# -*- coding: utf-8 -*-
# @Time     : 2020/4/26 10:20
# @Author   : zbc
# @Email    : zhangbc0315@outlook.com
# @FileName : rxn_reverse_predictor.py
# @Software : PyCharm

from typing import List, Dict, Iterator

import networkx as nx
from networkx.algorithms import isomorphism

from crkitpy.molgraph.molecule import Molecule
from crkitpy.molgraph.atom import Atom
from crkitpy.const import ATOM_PROP


class RxnReversePredictor:

    @classmethod
    def predict_by_rxn_template(cls, rxn_template: Molecule, product: Molecule) -> List[Molecule]:
        mapping_p_t_atoms = cls._match_template(rxn_template, product)
        if mapping_p_t_atoms == {}:
            return []
        cls._draw_product_mapping_result(product, mapping_p_t_atoms)
        cls._link_and_break_bonds(rxn_template, mapping_p_t_atoms, product)
        return list(product.get_detached_sub_mols())

    @classmethod
    def _draw_product_mapping_result(cls, product: Molecule, mapping_p_t_atoms: Dict[Atom, Atom]) -> None:
        for p_atom in mapping_p_t_atoms.keys():
            p_atom.set_property(ATOM_PROP.K_NEED_SAVE, True)
        product.draw()

    @classmethod
    def _link_and_break_bonds(cls, rxn_template, mapping_t_p_atoms, product):
        for used_broken_bond in rxn_template.get_used_broken_bonds():
            t_start_atom = used_broken_bond.get_atoms()[0]
            t_end_atom = used_broken_bond.get_atoms()[1]
            p_start_atom = cls._get_mapping_p_atom_by_t_atom(t_start_atom, mapping_t_p_atoms)
            p_end_atom = cls._get_mapping_p_atom_by_t_atom(t_end_atom, mapping_t_p_atoms)
            product.remove_bond(p_start_atom.get_bond(p_end_atom))
        for used_linked_bond in rxn_template.get_used_linked_bonds():
            p_atoms = []
            for t_atom in used_linked_bond.get_atoms():
                if t_atom.get_property(ATOM_PROP.K_MAYBE_RGROUP):
                    p_atoms.append(t_atom)
                    product.add_atom(t_atom)
                else:
                    p_atoms.append(cls._get_mapping_p_atom_by_t_atom(t_atom, mapping_t_p_atoms))
            product.add_bond_by_type(used_linked_bond.get_bond_type(), p_atoms[0], p_atoms[1])

    @staticmethod
    def _get_mapping_p_atom_by_t_atom(t_atom: Atom, mapping_p_t_atom: Dict[Atom, Atom]) -> Atom:
        for key, value in mapping_p_t_atom.items():
            if value == t_atom:
                return key

    @classmethod
    def _match_template(cls, product: Molecule, rxn_template: Molecule) -> Dict[Atom, Atom]:
        gm = isomorphism.GraphMatcher(rxn_template.get_mol_graph_for_match(), product.get_mol_graph_for_match(),
                                      node_match=cls._node_match_func,
                                      edge_match=cls._edge_match_func)
        res = gm.subgraph_is_isomorphic()
        return gm.mapping

    @staticmethod
    def _node_match_func(node, t_node) -> bool:
        if t_node['is_reacted']:
            return node['symbol'] == t_node['symbol'] and node['env'] == node['env']
        else:
            return node['symbol'] == t_node['symbol']

    @staticmethod
    def _edge_match_func(edge, t_edge) -> bool:
        return edge['bond_type'] == t_edge['bond_type']

