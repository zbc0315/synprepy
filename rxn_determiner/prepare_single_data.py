#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Time    : 2021/8/18 10:56
# @Author  : zhangbc0315@outlook.com
# @File    : prepare_single_data.py
# @Software: PyCharm

import sys
import os
from typing import Tuple, List, Iterator
import json

import pandas as pd
import torch

from nn_utils.gnn_data_utils import GnnDataUtils
from config import RDPath, SSPath


class PrepareSingleData:
    spectator_count_filter_fps = {'cat': SSPath.CATS_COUNT_FILTER_FP,
                                  'sol': SSPath.SOLS_COUNT_FILTER_FP}
    cat_count_filter_df = None
    sol_count_filter_df = None
    rxn_with_spectators_before_2016_df = None

    # ============================================
    # -- 过滤后的 催化剂/溶剂 统计数据 --

    @classmethod
    def get_spectators_by_rss(cls, spectator_type: str):
        """ 获得和 Reaction Spectators Selector 中一样的 Spectator Count Filter

        :param spectator_type:
        :return:
        """
        spectator_count_filter_fp = cls.spectator_count_filter_fps[spectator_type]
        df = pd.read_csv(spectator_count_filter_fp, sep='\t', encoding='utf-8')
        return df

    @classmethod
    def get_filter_spectators(cls, spectator_type: str):
        """ 获得 Spectator Count Filter，有两种途径：
        1. 使用和Reaction Spectators Selector项目中一样的数据，方便整体系统中的配合
        2. 从Reaction Spectator Count中重新过滤，在调试程序时方便查看不同过滤要求下的训练效果（未完成）

        :param spectator_type:
        :return:
        """
        return cls.get_spectators_by_rss(spectator_type)

    # ====================================================
    # -- 通过rid查询反应的催化剂/溶剂 --

    @classmethod
    def query_rxn_data_by_rid(cls, rid):
        idxes = cls.rxn_with_spectators_before_2016_df[
            (cls.rxn_with_spectators_before_2016_df.rid == rid)].index.tolist()
        if len(idxes) != 1:
            return None
        return cls.rxn_with_spectators_before_2016_df.loc[idxes[0], ]

    @classmethod
    def query_index_by_smiles(cls, df: pd.DataFrame, smiles: str):
        idxes = df[(df.smiles == smiles)].index.tolist()
        if len(idxes) != 1:
            return None
        return idxes[0]

    @classmethod
    def init_spectator_one_hot(cls, df):
        """ 初始化催化剂/溶剂的one-hot向量，假设一共有N个催化剂/溶剂，则one-hot向量长度为N+2
        0~N-1表示使用了已知的催化剂/溶剂
        N表示没有使用催化剂/溶剂
        N+1表示使用了未知的催化剂/溶剂

        :param df:
        :return:
        """
        res = [0] * (len(df) + 2)
        return res

    @classmethod
    def get_spectator_one_hot(cls, smileses: [str], df: pd.DataFrame) -> [int]:
        cats_num = len(df)
        res = cls.init_spectator_one_hot(df)
        if len(smileses) != 1:
            res[cats_num] = 1
            return res
        idx = cls.query_index_by_smiles(df, smileses[0])
        if idx is None:
            res[cats_num + 1] = 1
        else:
            res[idx] = 1
        return res

    @classmethod
    def get_spectators_one_hot_by_rid(cls, rid: int):
        rxn_data = cls.query_rxn_data_by_rid(rid)
        cats = eval(rxn_data['cats'])
        sols = eval(rxn_data['sols'])
        cat_one_hot = cls.get_spectator_one_hot(cats, cls.cat_count_filter_df)
        sol_one_hot = cls.get_spectator_one_hot(sols, cls.sol_count_filter_df)
        return torch.tensor(cat_one_hot, dtype=torch.float), \
               torch.tensor(sol_one_hot, dtype=torch.float)

    # ====================================================

    @classmethod
    def get_mol_gnn_data_from_inchi(cls, inchi: str, no_radical: bool = False):
        gnn_data = GnnDataUtils.gnn_data_from_inchi(inchi, no_radical=no_radical)
        return gnn_data

    @classmethod
    def get_mol_gnn_data_from_inchis(cls, inchis: [str]):
        gnn_data = GnnDataUtils.gnn_data_from_inchis(inchis)
        return gnn_data

    # ====================================================

    @classmethod
    def _get_wrong_rxn_from_one_file(cls, fp: str):
        """

        :param fp:
        :return: rid, wrong product inchi, reactants inchis
        """
        wr_df = pd.read_csv(fp, sep='\t', encoding='utf-8')
        for _, row in wr_df.iterrows():
            rid = row['rid']
            wp_inchi = row['wp_inchi']
            rs_inchis = json.loads(row['rs_inchis'])
            yield rid, wp_inchi, rs_inchis

    @classmethod
    def get_wrong_rxn(cls):
        """

        :return: rid, wrong product inchi, reactants inchis
        """
        num = 0
        for i in range(10):
            fp = RDPath.WRONG_RXN_FPS.format(i)
            print(fp)
            for rid, wp_inchi, rs_inchis in cls._get_wrong_rxn_from_one_file(fp):

                cat_ont_hot, sol_one_hot = cls.get_spectators_one_hot_by_rid(rid)
                wp_gnn_data = cls.get_mol_gnn_data_from_inchi(wp_inchi, no_radical=True)
                rs_gnn_data = cls.get_mol_gnn_data_from_inchis(rs_inchis)
                if wp_gnn_data is None or rs_gnn_data is None:
                    print(f"无法从 {rid} 得到gnn data")
                    continue
                data = {'cat': cat_ont_hot,
                        'sol': sol_one_hot,
                        'rs': rs_gnn_data,
                        'p': wp_gnn_data}
                gnn_fp = os.path.join(RDPath.RXN_WRONG_SINGLE_DATA_DP, f"w_{num}.pt")
                torch.save(data, gnn_fp)
                num += 1

    # ===========================================

    @classmethod
    def get_reaction_gnn_data_from_rid(cls, rid: int):
        rxn_data = cls.query_rxn_data_by_rid(rid)
        rs_smiles = json.loads(rxn_data['reactants'])
        ps_smiles = json.loads(rxn_data['products'])
        if len(rs_smiles) == 0 or len(ps_smiles) != 1:
            return None, None
        rs_data = GnnDataUtils.gnn_data_from_smileses(rs_smiles)
        p_data = GnnDataUtils.gnn_data_from_smiles(ps_smiles[0])
        return rs_data, p_data

    @classmethod
    def get_right_rxn(cls):
        rr_df = pd.read_csv(RDPath.RIGHT_RXN_FP, sep='\t', encoding='utf-8')
        num = 0
        error_num = 0
        for _, row in rr_df.iterrows():
            rid = row['rid']
            cat_smiles = row['cat']
            sol_smiles = row['sol']

            rs_gnn_data, p_gnn_data = cls.get_reaction_gnn_data_from_rid(rid)
            if rs_gnn_data is None or p_gnn_data is None:
                error_num += 1
                print(f"{num} - {error_num} 无法从{rid}获得gnn data")
                continue

            cat_one_hot, sol_one_hot = cls.get_spectators_one_hot_by_rid(rid)
            # sol_one_hot = cls.get_spectators_one_hot_by_rid(rid)

            res = {'cat': cat_one_hot,
                   'sol': sol_one_hot,
                   'rs': rs_gnn_data,
                   'p': p_gnn_data}
            gnn_fp = os.path.join(RDPath.RXN_RIGHT_SINGLE_DATA_DP, f"r_{num}.pt")
            torch.save(res, gnn_fp)
            num += 1

    @classmethod
    def init(cls):
        cls.cat_count_filter_df = cls.get_filter_spectators('cat')
        cls.sol_count_filter_df = cls.get_filter_spectators('sol')
        cls.rxn_with_spectators_before_2016_df = pd.read_csv(SSPath.RXN_WITH_SPECTATORS_FP, sep='\t', encoding='utf-8')

    @classmethod
    def process(cls, rxn_type: str):
        cls.init()
        if rxn_type == 'right':
            cls.get_right_rxn()
        else:
            cls.get_wrong_rxn()

    @classmethod
    def test(cls):
        cls.init()

        res = cls.query_rxn_data_by_rid(0)
        print(res)
        res = cls.query_rxn_data_by_rid(2140633)
        print(res)

        res = cls.get_spectators_one_hot_by_rid(2140633)
        print(res)


if __name__ == "__main__":
    rxn_type = sys.argv[1]
    PrepareSingleData.process(rxn_type)
