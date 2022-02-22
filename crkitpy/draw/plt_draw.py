# -*- coding: utf-8 -*-
# @Time     : 2020/4/21 10:46
# @Author   : zbc
# @Email    : zhangbc0315@outlook.com
# @FileName : plt_draw.py
# @Software : PyCharm

from typing import Dict, List

import networkx as nx
import matplotlib.pyplot as plt


class PltDraw:

    @staticmethod
    def draw(nx_graph: nx.Graph, width, node_color) -> None:
        positions = nx.kamada_kawai_layout(nx_graph)
        # positions = nx.fruchterman_reingold_layout(nx_graph)
        nx.draw(nx_graph,
                pos=positions,
                node_size=300,
                node_color=node_color,
                font_color='black',
                with_labels=True)
        nx.draw_networkx_edges(nx_graph, pos=positions, width=width)

    @staticmethod
    def draw_graphs(graph_dict: Dict[str, List[nx.Graph]]) -> None:
        row_count = len(graph_dict.keys())
        col_count = max([len(graph_dict[key]) for key in graph_dict.keys()])
        for row_num, key in enumerate(graph_dict.keys()):
            for col_num, nx_graph in enumerate(graph_dict[key]):
                plot_num = row_num*row_count + col_num + 1
                plt.subplot(row_count, col_count, plot_num, label=key)
                nx.draw(nx_graph,
                        pos=nx.kamada_kawai_layout(nx_graph),
                        node_size=300,
                        node_color='w',
                        font_color='black',
                        with_labels=True)
        plt.show()

    @staticmethod
    def show():
        plt.show()
