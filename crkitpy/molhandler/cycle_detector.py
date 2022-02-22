# -*- coding: utf-8 -*-
# @Time     : 2020/5/27 11:55
# @Author   : zbc
# @Email    : zhangbc0315@outlook.com
# @FileName : cycle_detector.py
# @Software : PyCharm

from typing import List, Set

import matplotlib.pyplot as plt
import networkx as nx

from crkitpy.molgraph.molecule import Molecule


class CycleDetector:

    # ================
    # -- Get Cycles --

    @classmethod
    def get_cycles(cls, molecule: Molecule, is_right_node=None):
        cycles = []
        passed_nodes = set()
        graph = molecule.mol_graph
        branch_nodes = cls._get_branch_nodes(graph)
        for (node1, node2) in graph.edges:
            if node1 in branch_nodes or node2 in branch_nodes:
                continue
            if node1 in passed_nodes and node2 in passed_nodes:
                continue
            if is_right_node is not None:
                if not (is_right_node(node1) and is_right_node(node2)):
                    continue
            passed_nodes.add(node1)
            passed_nodes.add(node2)
            edge_data = graph.get_edge_data(node1, node2)
            graph.remove_edge(node1, node2)
            try:
                cycle = nx.dijkstra_path(graph, node1, node2)
                if is_right_node is None or all([is_right_node(n) for n in cycle]):
                    cycles.append(cycle)
                    passed_nodes = passed_nodes.union(cycle)
            except Exception as e:
                pass
            graph.add_edge(node1, node2, bond_type=edge_data['bond_type'], bond_to_r=['bond_to_r'],
                                bond_weight=edge_data['bond_weight'])
        return cycles

    @classmethod
    def _get_branch_nodes(cls, graph: nx.Graph) -> Set:
        """ 获得肯定不在环上的节点序列

        :param graph: 图
        :return: 肯定不在环上的节点序列
        """
        res = set()
        res_len = len(res)
        while True:
            for node in graph.nodes:
                neighbors = nx.neighbors(graph, node)
                neighbors_num = len(set(neighbors) - res)
                if neighbors_num <= 1:
                    res.add(node)
            if res_len == len(res):
                break
            res_len = len(res)
        return res

    @classmethod
    def _get_cycles(cls, graph: nx.Graph):
        tree = nx.Graph()
        nodes_parents = {}
        first_node = list(graph.nodes)[0]
        tree.add_node(first_node)
        cycles = []
        nodes = [first_node]
        while len(tree.nodes) < len(graph.nodes):
            new_nodes = []
            for node in nodes:
                child_nodes = cls._get_child_nodes(graph, node, cls._get_parent_node(nodes_parents, node))
                for child_node in child_nodes:
                    if tree.has_node(child_node):
                        cycles.append(nx.dijkstra_path(tree, node, child_node))
                        continue
                    else:
                        tree.add_node(child_node)
                        tree.add_edge(node, child_node)
                        nodes_parents[child_node] = node
                        new_nodes.append(child_node)
            # nx.draw(tree, with_labels=True)
            # plt.show()
            nodes = new_nodes

        # cls._parse_node(first_node, graph, tree, nodes_parents, cycles)
        return cycles

    @classmethod
    def _parse_node(cls, node, graph, tree, nodes_parents, cycles):
        child_nodes = cls._get_child_nodes(graph, node, cls._get_parent_node(nodes_parents, node))
        cycled_nodes = []
        for child_node in child_nodes:
            if tree.has_node(child_node):
                cycles.append(nx.dijkstra_path(tree, node, child_node))
                cycled_nodes.append(child_node)
                continue
            else:
                tree.add_node(child_node)
                tree.add_edge(node, child_node)
                nodes_parents[child_node] = node
                # nx.draw(tree, with_labels=True)
                # plt.show()
        if len(tree.nodes) == len(graph.nodes):
            return
        else:
            for child_node in child_nodes:
                if child_node not in cycled_nodes:
                    cls._parse_node(child_node, graph, tree, nodes_parents, cycles)

    @classmethod
    def _get_parent_node(cls, nodes_parents: dict, node):
        if node in nodes_parents:
            return nodes_parents[node]
        else:
            return None

    @classmethod
    def _get_child_nodes(cls, graph: nx.Graph, node, parent_node):
        neighbor_nodes = list(graph.neighbors(node))
        if parent_node is not None:
            neighbor_nodes.remove(parent_node)
        return neighbor_nodes
