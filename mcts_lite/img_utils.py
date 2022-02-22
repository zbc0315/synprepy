#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Time    : 2021/10/17 22:35
# @Author  : zhangbc0315@outlook.com
# @File    : img_utils.py
# @Software: PyCharm

import os

import matplotlib.pyplot as plt
import pandas as pd
import networkx as nx
from networkx.drawing.nx_agraph import graphviz_layout


class ImgUtils:

    _unsolved = '#193498'
    _solved = '#4E9F3D'

    _edge_color = []
    _node_color = []
    _node_labels = [0]

    # _nid = 11586
    _nid = 8846
    _debug_dp = '/mnt/c/Users/zhang/Documents/Data/RxnMolPredictor/Mcts/debug/'
    _debug_dp = 'C:\\Users\\zhang\\Documents\\Data\\RxnMolPredictor\\Mcts\\debug'
    _rxn_fp = os.path.join(_debug_dp, f'{_nid}-Rxn.tsv')
    _rxn_df = pd.read_csv(_rxn_fp, sep='\t', encoding='utf-8')
    _node_fp = os.path.join(_debug_dp, f'{_nid}-Node.tsv')
    _node_df = pd.read_csv(_node_fp, sep='\t', encoding='utf-8')

    @classmethod
    def set_node_edge_color(cls, tree, solved_nids):
        nodes = tree.nodes()
        edges = tree.edges()
        cls._node_color = [cls._solved if nid in solved_nids else cls._unsolved for nid in nodes]
        cls._edge_color = [cls._solved if edge[0] in solved_nids and edge[1] in solved_nids else cls._unsolved for edge in edges]

    @classmethod
    def get_tid_from_nid(cls, nid):
        tid = int(cls._rxn_df.tid[cls._node_df.rid[nid]] / 100)
        return f'{tid}'

    @classmethod
    def get_solved_nodes(cls):
        for current_nid in cls._node_df.query(f'is_solved').index:
            while True:
                yield int(current_nid)
                if current_nid == 0:
                    break
                current_nid = cls._node_df.parent_nid[current_nid]

    @classmethod
    def create_mct(cls, tree, current_nid):
        for n, child_nid in enumerate(eval(cls._node_df.child_nids[current_nid])):
            child_nid = int(child_nid)
            # if child_nid in [2, 3, 4, 5, 12, 13, 14, 20, 21, 23, 24, 26, 27, 28, 29, 35, 37, 38, 42, 48, 70, 102, 926, 927, 932, 936, 1109, 1110, 1112, 3503, 3959, 5295, 6101, 6444, 6, 7, 8, 9, 10, 11]:
            #     continue
            # if 17 >= child_nid >= 5:
            #     continue
            if len(eval(cls._node_df.child_nids[child_nid])) == 0 and not cls._node_df.is_solved[child_nid]:
                continue
            # if n > 10:
            #     break
            tid = cls.get_tid_from_nid(child_nid)
            # tree.add_node(child_nid, desc=f'{child_nid}')
            tree.add_node(child_nid, desc=f'{tid}')
            tree.add_edge(current_nid, child_nid)
            cls._node_labels.append(n)
            cls.create_mct(tree, child_nid)
        return tree

    @classmethod
    def draw_mct(cls):
        solved_nids = list(cls.get_solved_nodes())
        tree = nx.DiGraph()
        tree.add_node(0, desc='root')
        edge_color = []
        tree = cls.create_mct(tree, 0)
        pos = graphviz_layout(tree, prog='dot')
        # pos = nx.spiral_layout(tree)
        # plt.figure(figsize=(10, 10), dpi=200)
        plt.figure(figsize=(10, 10))
        cls.set_node_edge_color(tree, solved_nids)
        # nx.draw_networkx_nodes(tree, pos, nodelist=tree.nodes(), node_color=cls._node_color, node_shape='o', node_size=800)
        # nx.draw_networkx_edges(tree, pos, edgelist=tree.edges())
        nx.draw(tree, pos, node_size=1500, arrows=False, width=5, node_shape='o', node_color=cls._node_color, edge_color=cls._edge_color, with_labels=False)
        # nx.draw(tree, pos, node_size=0, arrows=False, width=5, node_shape='o', node_color=cls._node_color, edge_color=cls._edge_color, with_labels=False)
        # node_labels = nx.get_node_attributes(tree, 'desc')
        # nx.draw_networkx_labels(tree, pos, labels=node_labels, font_color='w', font_size=11)
        plt.axis('equal')
        plt.show()

    @classmethod
    def draw_mct_from_fp(cls, node_fp):
        cls.draw_mct()


if __name__ == "__main__":
    # fp = 'C:\\Users\\zhang\\Documents\\Data\\RxnMolPredictor\\Mcts\\debug\\2271-Node.tsv'
    fp = '/mnt/c/Users/zhang/Documents/Data/RxnMolPredictor/Mcts/debug/1637-Node.tsv'
    ImgUtils.draw_mct_from_fp(fp)
