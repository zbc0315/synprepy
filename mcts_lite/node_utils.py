#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Time    : 2021/9/13 13:39
# @Author  : zhangbc0315@outlook.com
# @File    : node_utils.py
# @Software: PyCharm
import math

import pandas as pd

from mcts_lite.mol_utils import MolUtils
from mcts_lite.rxn_utils import RxnUtils


class NodeUtils:

    node_df = None
    max_level = 25
    c = 0.5

    @classmethod
    def save(cls, fp: str):
        if cls.node_df is not None:
            cls.node_df.to_csv(fp, sep='\t', encoding='utf-8')

    @classmethod
    def init(cls, max_level):
        cls.node_df = pd.DataFrame(columns=['level', 'num', 'mids', 'rid',
                                            'is_solved', 'no_expand', 'no_select', 'parent_nid', 'child_nids',
                                            'reward', 'Q', 'P', 'visited_times'])
        cls.max_level = max_level

    @classmethod
    def get_rest_step(cls, nid):
        level = int(cls.node_df.level[nid])
        return cls.max_level - level

    # ==========
    # -- root --

    @classmethod
    def create_root_node_by_inchi(cls, target_inchi):
        target_mid = MolUtils.add_or_query_inchi(target_inchi)
        nid = cls.create_root_node_by_mid(target_mid)
        return nid

    @classmethod
    def create_root_node_by_mid(cls, target_mid):
        cls._insert_node(0, 0, {target_mid}, None, [], None, P=1)
        return 0

    # =====================
    # -- parent children --

    @classmethod
    def get_parent_node_data(cls, node_data):
        return cls.query_by_nid(node_data.parent_nid)

    @classmethod
    def get_child_nodes_data(cls, node_data):
        for child_nid in node_data.child_nids:
            yield cls.query_by_nid(child_nid)

    # =======================
    # -- create child node --

    @classmethod
    def _update_mids(cls, mids, rid):
        res = mids.copy()
        rxn_data = RxnUtils.query_by_rid(rid)
        res = res - {rxn_data.p_mid}
        res = res.union(rxn_data.rs_mids)
        return res

    @classmethod
    def _get_is_solved(cls, mids):
        return all([MolUtils.query_by_mid(mid).is_cac for mid in mids])

    @classmethod
    def create_child_nodes_by_rids(cls, nid, rids):
        node_data = cls.query_by_nid(nid)
        for n, rid in enumerate(rids):
            P = cls.get_P_by_rid(rid)
            mids = cls._update_mids(node_data.mids, rid)
            is_solved = cls._get_is_solved(mids)
            child_nid = cls._insert_node(level=node_data.level+1,
                                         num=n,
                                         mids=mids,
                                         rid=rid,
                                         parent_nid=nid,
                                         child_nids=[],
                                         P=math.pow(P, 1),
                                         is_solved=is_solved)
            node_data.child_nids.append(child_nid)

    # ===========
    # -- query --

    @classmethod
    def query_by_nid(cls, nid):
        return cls.node_df.loc[nid, ]

    @classmethod
    def is_leaf(cls, nid):
        return len(cls.query_by_nid(nid).child_nids) == 0

    @classmethod
    def is_terminal(cls, nid):
        node_data = cls.query_by_nid(nid)
        return node_data.is_solved or node_data.no_expand or node_data.no_select

    @classmethod
    def get_un_terminal_child_nids(cls, nid):
        node_data = NodeUtils.query_by_nid(nid)
        for child_nid in node_data.child_nids:
            if not cls.is_terminal(child_nid):
                yield child_nid

    @classmethod
    def get_unsolved_inchis_and_mids(cls, nid):
        node_data = cls.query_by_nid(nid)
        for mid in node_data.mids:
            mol_data = MolUtils.query_by_mid(mid)
            if not mol_data.is_cac:
                yield mol_data.inchi, mid

    @classmethod
    def get_unsolved_inchi_and_mid(cls, nid):
        inchis_mids = list(cls.get_unsolved_inchis_and_mids(nid))
        inchis_mids_sorted = sorted(inchis_mids, key=lambda x: len(x[0]))
        return inchis_mids_sorted[-1]

    # =======================
    # -- get past products --

    @classmethod
    def _get_intermediate_product(cls, nid, mid):
        """ 获得mid为原料时，产生的产物

        """
        node = cls.query_by_nid(nid)
        nid = nid
        while node.level > 0:
            rxn = RxnUtils.query_by_rid(node.rid)
            if mid in rxn.rs_mids:
                return rxn.p_mid, nid
            nid = node.parent_nid
            node = cls.query_by_nid(nid)
        return None, None

    @classmethod
    def get_intermediate_products(cls, nid, mid):
        """ 获得从mid到目标产物的合成路径中的所有中间产物，
        要求expand阶段，mid的原料中不能再出现这些中间产物

        """
        res = []
        nid = nid
        mid = mid
        while mid is None:
            mid, nid = cls._get_intermediate_product(nid, mid)
            # res.append(MolUtils.query_by_mid(mid).inchi)
            res.append(mid)
        return res

    # ===========================

    # @classmethod
    # def get_past_products(cls, nid):
    #     res = []
    #     node = cls.query_by_nid(nid)
    #     while node.level > 0:
    #         rxn = RxnUtils.query_by_rid(node.rid)
    #         p_mol = MolUtils.query_by_mid(rxn.p_mid)
    #         res.append(p_mol.inchi)
    #
    #         node = NodeUtils.query_by_nid(node.parent_nid)
    #     return res

    # =====================
    # -- calculate score --

    @classmethod
    def get_P_by_rid(cls, rid):
        rxn_data = RxnUtils.query_by_rid(rid)
        return rxn_data.prob

    @classmethod
    def update_Q(cls, nid):
        node_data = cls.query_by_nid(nid)

        XI = node_data.level
        parent_node_data = cls.get_parent_node_data(node_data)
        while parent_node_data.level > 0:
            XI -= parent_node_data.P * 0.99
            parent_node_data = cls.get_parent_node_data(parent_node_data)

        W = max(0.0, (cls.max_level - XI)/cls.max_level)

        Q = 0
        for child_node_data in cls.get_child_nodes_data(node_data):
            Q += child_node_data.reward * W
        Q = Q / node_data.visited_times

        cls.change_Q(Q, nid)

    @classmethod
    def get_A(cls, nid):
        node_data = cls.query_by_nid(nid)
        N = node_data.visited_times
        N_1 = cls.get_parent_node_data(node_data).visited_times
        return node_data.Q / N + cls.c * node_data.P * math.sqrt(N_1) / (1 + N)

    # =====================
    # -- insert / change --

    @classmethod
    def change_Q(cls, Q, nid):
        cls.node_df.loc[nid, 'Q'] = Q

    @classmethod
    def visit(cls, nid):
        cls.node_df.loc[nid, 'visited_times'] += 1

    @classmethod
    def set_is_solved(cls, nid, is_solved):
        cls.node_df.loc[nid, 'is_solved'] = is_solved

    @classmethod
    def set_no_expand(cls, nid, no_expand):
        cls.node_df.loc[nid, 'no_expand'] = no_expand

    @classmethod
    def set_no_select(cls, nid, no_select):
        cls.node_df.loc[nid, 'no_select'] = no_select

    @classmethod
    def set_reward(cls, nid, reward):
        cls.node_df.loc[nid, 'reward'] = reward

    @classmethod
    def _insert_node(cls, level, num, mids, parent_nid, child_nids, rid,
                     is_solved=False, no_expand=False, no_select=False,
                     reward=0, Q=0, P=None, visited_times=1):
        cls.node_df = cls.node_df.append({'level': level, 'num': num, 'mids': mids, 'rid': rid,
                                          'parent_nid': parent_nid, 'child_nids': child_nids,
                                          'is_solved': is_solved, 'no_expand': no_expand, 'no_select': no_select,
                                          'reward': reward, 'Q': Q, 'P': P, 'visited_times': visited_times},
                                         ignore_index=True)
        return len(cls.node_df) - 1

    @classmethod
    def get_rxn_plan(cls, nid):
        """
        化学反应的smiles列表
        化学反应的inchi字典列表
        用到的原料
        """
        res_smiles = []
        res_dict = []
        current_nid = nid
        current_node = NodeUtils.query_by_nid(current_nid)
        building_blocks = [MolUtils.query_by_mid(mid).inchi for mid in current_node.mids]
        while current_node.level > 0:
            res_smiles.append(RxnUtils.get_rxn_smiles(current_node.rid))
            res_dict.append(RxnUtils.get_rxn_dict(current_node.rid))
            current_nid = current_node.parent_nid
            current_node = NodeUtils.query_by_nid(current_nid)
        return res_smiles, res_dict, building_blocks


if __name__ == "__main__":
    NodeUtils.init(25)
    NodeUtils.create_root_node_by_mid(0)
    r = NodeUtils.query_by_nid(0)
