# -*- coding: utf-8 -*-
# @Time     : 2020/4/21 16:56
# @Author   : zbc
# @Email    : zhangbc0315@outlook.com
# @FileName : rxn_forward_predictor.py
# @Software : PyCharm

from typing import List, Dict

import networkx as nx
from networkx.algorithms import isomorphism

from crkitpy.molgraph.reaction import Reaction
from crkitpy.molgraph.molecule import Molecule
from crkitpy.molgraph.atom import Atom
from crkitpy.utils.array_utils import ArrayUtils
from crkitpy.const import BOND_TYPE, ATOM_PROP


class RxnForwardPredictor:

    @classmethod
    def predict_by_rxn_template(cls, rxn_template: Reaction, reactants: List[Molecule]) -> List[List[Molecule]]:
        """ 使用化学反应模板预测化学反应

        :param rxn_template: 化学反应模板
        :param reactants: 待预测产物的反应物
        :return: 产物列表的列表
        """
        rs = []
        for reactant in reactants:
            r, _, __ = reactant.copy()
            r.clear_mapping()
            rs.append(r)
        # TODO 理应有多组matched_reactants，并产生多个可能的产物
        matched_reactants = cls._match_reactants(rxn_template, rs)
        cls._map_reactants(matched_reactants)
        cls._remove_template_atoms_in_reactant(matched_reactants)
        products = []
        for p_template in rxn_template.get_products():
            products.append(cls._get_product(p_template, rs))
        return [products]

    @classmethod
    def _get_product(cls, p_template: Molecule, reactants: List[Molecule]) -> Molecule:
        """ 根据产物模板得到产物

        :param p_template: 产物模板
        :param reactants: 带预测产物的反应物
        :return: 预测的产物结果
        """
        for reactant in reactants:
            p_template.compose_graph(reactant)
            for r_link_atom in reactant.get_atoms_by_map_num(-1):
                if r_link_atom.get_property(ATOM_PROP.K_NEED_REMOVE) == 0:
                    continue
                p_link_atoms = list(p_template.get_atoms_by_map_num(r_link_atom.get_property(ATOM_PROP.K_NEED_REMOVE)))
                if len(p_link_atoms) != 1:
                    raise ValueError("Error Product Template With no Map Atom or Duplicated Map Num")
                p_template.add_bond_by_type(r_link_atom.get_property(ATOM_PROP.K_BROCK_BOND),
                                            r_link_atom,
                                            p_link_atoms[0])
        all_atoms = p_template.get_atoms().copy()
        for p_atom in all_atoms:
            if p_atom.get_symbol() == '*':
                p_template.remove_atom(p_atom)
            elif len(p_atom.get_bonds()) == 0:
                p_template.remove_atom(p_atom)
        return p_template

    @staticmethod
    def _remove_template_atoms_in_reactant(matched_reactants: Dict[Molecule, Molecule]):
        """ 移除反应物中与反应物模板中原子对应的原子

        :param matched_reactants: 反应物模板与反应物的匹配字典
        :return: None
        """
        for r_template, reactant in matched_reactants.items():
            for atom in r_template.get_atoms():
                if atom.get_symbol == '*':
                    continue
                if atom.map_num is None:
                    raise ValueError('Error Atom in Reactant Template with None Map Num')
                for need_removed_atom in reactant.get_atoms_by_map_num(atom.map_num):
                    reactant.remove_atom(need_removed_atom)
            pass

    @staticmethod
    def _get_r_group_neighbor(r_group: Atom) -> Atom:
        """ 得到'*'原子的临近原子，当前默认'*'原子仅与一个原子临近

        :param r_group: '*'原子
        :return: '*'原子的临近原子
        """
        neighbors = list(r_group.get_neighbors())
        if len(neighbors) != 1:
            raise ValueError("Error * atom with {} neighbors".format(len(neighbors)))
        return neighbors[0]

    @classmethod
    def _map_reactants(cls, matched_reactants: Dict[Molecule, Molecule]) -> None:
        """ 以反应物匹配到的反应物模板为基准，对反应物中的原子赋值map_num，map_num的值与匹配的反应物模板中对应的原子的map_num相同

        :param matched_reactants: 反应物模板与反应物的匹配字典
        :return: None
        """
        for r_template, reactant in matched_reactants.items():
            # TODO ERROR r_template 与 reactant 可能弄反了
            atoms_map = cls._get_map_of_match_rxn(r_template, reactant)
            for reactant_atom, template_atom in atoms_map.items():
                if template_atom.get_symbol() == '*':
                    need_remove = cls._get_r_group_neighbor(template_atom)
                    if need_remove.map_num == 15:
                        print(need_remove)
                    reactant_atom.set_property(ATOM_PROP.K_NEED_REMOVE, need_remove.map_num)
                    reactant_atom.set_property(ATOM_PROP.K_BROCK_BOND,
                                               r_template.get_bond(reactant_atom,
                                                                   cls._get_atom_from_atoms_map_by_template_atom(atoms_map, need_remove)).get_bond_type())
                    reactant_atom.map_num = -1
                else:
                    reactant_atom.map_num = template_atom.map_num
            for r_atom in reactant.get_atoms():
                for neighbor_atom in r_atom.get_neighbors():
                    if r_atom.map_num is not None or neighbor_atom.map_num is None or neighbor_atom.map_num == -1:
                        continue
                    i = neighbor_atom.map_num
                    if neighbor_atom.map_num == 15:
                        print(i)
                    r_atom.map_num = -1
                    r_atom.set_property(ATOM_PROP.K_NEED_REMOVE, neighbor_atom.map_num)
                    r_atom.set_property(ATOM_PROP.K_BROCK_BOND, reactant.get_bond(r_atom, neighbor_atom).get_bond_type())

    @staticmethod
    def _get_atom_from_atoms_map_by_template_atom(atoms_map: Dict[Atom, Atom], template_atom: Atom) -> Atom:
        for r_atom, t_atom in atoms_map.items():
            if template_atom == t_atom:
                return r_atom
        raise ValueError('cannot get atom from atoms map by template atom')

    @classmethod
    def _match_reactants(cls, rxn_template: Reaction, reactants: List[Molecule]) -> Dict[Molecule, Molecule]:
        """ 匹配反应物模板与反应物

        :param rxn_template: 化学反应模板，里面包含一个或多个反应物模板及一个或多个产物模板
        :param reactants: 待预测产物的反应物
        :return: 匹配的反应物模板与产物模板的字典
        """
        temp_matched_reactants = {}
        for r_template in rxn_template.get_reactants():
            temp_matched_reactants[r_template] = []
            for reactant in reactants:
                if cls._is_matched_rxn(r_template, reactant):
                    temp_matched_reactants[r_template].append(reactant)
        matched_reactants = ArrayUtils.select_unique_array_in_matrix(temp_matched_reactants)
        if matched_reactants == {}:
            matched_reactants = ArrayUtils.select_random_array_in_matrix(temp_matched_reactants)
        return matched_reactants

    @classmethod
    def _is_matched_rxn(cls, r_template: Molecule, reactant: Molecule) -> bool:
        """ 判断反应物模板是否与反应物的子图同构

        :param r_template: 反应物模板
        :param reactant: 反应物
        :return: bool（是否匹配）
        """
        gm = cls._match_rxn(r_template, reactant)
        return gm.subgraph_is_isomorphic()

    @classmethod
    def _get_map_of_match_rxn(cls, r_template: Molecule, reactant: Molecule) -> Dict[Atom, Atom]:
        """ 得到反应物模板与反应物中一一对应的原子的对应关系
        此处有点多余，在判断反应物模板是否与反应物的子图同构时就应该得到原子的对应关系，
        如有必要，可通过修改utils.array_utils.py及本文件中部分代码省去本步骤

        :param r_template: 反应物模板
        :param reactant: 反应物
        :return: 反应物模板与反应物的原子的对应关系，其中，key是反应物中的原子，value是反应物模板中的原子
        """
        gm = cls._match_rxn(r_template, reactant)
        return gm.mapping

    @classmethod
    def _match_rxn(cls, r_template: Molecule, reactant: Molecule) -> isomorphism.GraphMatcher:
        """ 匹配反应物模板与反应物

        :param r_template: 反应物模板
        :param reactant: 反应物
        :return: 匹配结果，networkx.algorithms.isomorphism.GraphMatcher
        """
        gm = isomorphism.GraphMatcher(reactant.get_mol_graph_for_match(), r_template.get_mol_graph_for_match(),
                                      node_match=cls._node_match_func,
                                      edge_match=cls._edge_match_func)
        gm.subgraph_is_isomorphic()
        return gm

    @staticmethod
    def _node_match_func(node, template_node):
        """ 匹配反应物模板与反应物时的节点匹配原则
            节点匹配需要考虑节点本身的symbol以及节点周围节点的symbol
            例外：'*'节点与任意节点匹配

        :param node: 反应物中的节点
        :param template_node: 反应物模板中的节点
        :return: bool（是否匹配）
        """
        if node['symbol'] == '*' or template_node['symbol'] == '*':
            return True
        return node['symbol'] == template_node['symbol'] and node['env'] == node['env']

    @staticmethod
    def _edge_match_func(edge, template_edge):
        """ 匹配反应物模板与反应物时的边缘匹配原则
            边缘匹配需要考虑边缘的bond_type
            例外：反应物模板中与'*'节点相连的边缘的bond_type如果是SINGLE或DOUBLE，在对应的反应物模板中允许匹配到bond_type为AROMATIC的边缘

        :param edge: 反应物中的边缘
        :param template_edge: 反应物模板中的边缘
        :return: bool（是否匹配）
        """

        if template_edge['bond_to_r'] and template_edge['bond_type'] in [BOND_TYPE.SINGLE, BOND_TYPE.DOUBLE] and edge['bond_type'] == BOND_TYPE.AROMATIC:
            return True
        return edge['bond_type'] == template_edge['bond_type']
