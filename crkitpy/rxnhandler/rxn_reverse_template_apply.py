# -*- coding: utf-8 -*-
# @Time     : 2020/5/29 9:40
# @Author   : zbc
# @Email    : zhangbc0315@outlook.com
# @FileName : rxn_reverse_template_apply.py
# @Software : PyCharm

from typing import List, Dict, Iterator

import networkx as nx
from networkx.algorithms import isomorphism

from crkitpy.molgraph.molecule import Molecule
from crkitpy.molgraph.atom import Atom
from crkitpy.const import ATOM_PROP
from crkitpy.const import BOND_TYPE
from crkitpy.errors.bond_error import *


class RxnReverseTemplateApply:

    @classmethod
    def predict_by_rxn_template(cls, rxn_template: Molecule, product: Molecule) -> Iterator[List[Molecule]]:
        for mapping_t_p_atoms in cls._match_template(rxn_template, product):
            try:
                copy_product, p_old_new_atoms, _ = product.copy()
                copy_template, t_old_new_atoms, _ = rxn_template.copy()
                mapping_t_p_atoms = cls._reset_mapping_t_p_atoms(mapping_t_p_atoms, p_old_new_atoms, t_old_new_atoms)
                yield cls._reverse_react(copy_template, copy_product, mapping_t_p_atoms)
            except Exception as e:
                raise e

    @classmethod
    def _reset_mapping_t_p_atoms(cls, mapping_t_p_atoms: Dict[Atom, Atom],
                                 p_old_new_atoms: Dict[Atom, Atom],
                                 t_old_new_atoms: Dict[Atom, Atom]):
        new_mapping_t_p_atoms = {}
        for t_atom, p_atom in mapping_t_p_atoms.items():
            new_mapping_t_p_atoms[t_old_new_atoms[t_atom]] = p_old_new_atoms[p_atom]
        return new_mapping_t_p_atoms

    @classmethod
    def _reverse_react(cls, rxn_template: Molecule, product: Molecule, mapping_t_p_atoms: Dict[Atom, Atom]):
        cls._add_invisible_atoms(rxn_template, product)
        if mapping_t_p_atoms == {}:
            return []
        # cls._draw_product_mapping_result(product, mapping_t_p_atoms)
        cls._link_and_break_bonds(rxn_template, mapping_t_p_atoms, product)
        cls._change_charge(rxn_template, mapping_t_p_atoms)
        cls._changed_valence(rxn_template, mapping_t_p_atoms)
        cls._correct_aromatic(product)
        return list(product.get_detached_sub_mols())

    @classmethod
    def _correct_aromatic(cls, product: Molecule) -> None:
        for bond in product.get_bonds():
            if bond.get_bond_type() == BOND_TYPE.AROMATIC:
                for atom in bond.get_atoms():
                    atom.is_aromatic = True

    @classmethod
    def _change_charge(cls, rxn_template: Molecule, mapping_t_p_atoms: Dict[Atom, Atom]):
        for t_atom, (old_charge, new_charge) in rxn_template.get_changed_charge().items():
            p_atom = mapping_t_p_atoms[t_atom]
            p_atom.charge = old_charge

    @classmethod
    def _changed_valence(cls, rxn_template: Molecule, mapping_t_p_atoms: Dict[Atom, Atom]):
        for t_atom, (old_valence, new_valence) in rxn_template.get_changed_valence().items():
            p_atom = mapping_t_p_atoms[t_atom]
            p_atom.set_implicit_valence(old_valence - p_atom.get_explicit_valence())

    @classmethod
    def _add_invisible_atoms(cls, rxn_template: Molecule, product: Molecule):
        for atom in rxn_template.get_atoms():
            if not atom.is_visible:
                product.add_atom(atom)
        for bond in rxn_template.get_bonds():
            if all([not a.is_visible for a in bond.get_atoms()]):
                product.add_bond_by_type(bond.get_bond_type(), bond.get_atoms()[0], bond.get_atoms()[1])

    # ==========================
    # -- link and break bonds --

    @classmethod
    def _link_and_break_bonds(cls, rxn_template, mapping_t_p_atoms, product):
        for used_broken_bond in rxn_template.get_used_broken_bonds():
            t_start_atom = used_broken_bond.get_atoms()[0]
            t_end_atom = used_broken_bond.get_atoms()[1]
            p_start_atom = mapping_t_p_atoms[t_start_atom]
            p_end_atom = mapping_t_p_atoms[t_end_atom]
            product.remove_bond(p_start_atom.get_bond(p_end_atom))
        for used_linked_bond in rxn_template.get_used_linked_bonds():
            p_atoms = []
            for t_atom in used_linked_bond.get_atoms():
                if not t_atom.is_visible:
                    p_atoms.append(t_atom)
                else:
                    p_atoms.append(mapping_t_p_atoms[t_atom])
            product.add_bond_by_type(used_linked_bond.get_bond_type(), p_atoms[0], p_atoms[1])

    # ====================
    # -- match template --

    @classmethod
    def _match_template(cls, rxn_template: Molecule, product: Molecule) -> Iterator[Dict[Atom, Atom]]:
        gm = isomorphism.GraphMatcher(product.get_mol_graph_visible(),
                                      rxn_template.get_mol_graph_visible(),
                                      node_match=cls._node_match_func,
                                      edge_match=cls._edge_match_func)
        res = gm.subgraph_isomorphisms_iter()
        for mapping in res:
            yield dict(zip(mapping.values(), mapping.keys()))
#        return dict(zip(gm.mapping.values(), gm.mapping.keys()))

    @classmethod
    def match_sub_graph(cls, sub_graph: nx.Graph, graph: nx.Graph) -> Dict[Atom, Atom]:
        gm = isomorphism.GraphMatcher(graph, sub_graph,
                                      node_match=cls._node_match_func,
                                      edge_match=cls._edge_match_func)
        res = gm.subgraph_is_isomorphic()
        return gm.mapping

    @staticmethod
    def _node_match_func(node, t_node) -> bool:
        return node['symbol'] == t_node['symbol'] and \
               node['charge'] == t_node['charge'] and \
               node['is_aromatic'] == t_node['is_aromatic'] and \
               node['valence'] == t_node['valence']  # TODO 解析 template 时难以判断valence错误

    @staticmethod
    def _edge_match_func(edge, t_edge) -> bool:
        return edge['bond_type'] == t_edge['bond_type'] \
               or (edge['bond_type'] == BOND_TYPE.AROMATIC and t_edge['bond_type'] == BOND_TYPE.DOUBLE) \
               or (edge['bond_type'] == BOND_TYPE.DOUBLE and t_edge['bond_type'] == BOND_TYPE.AROMATIC)
