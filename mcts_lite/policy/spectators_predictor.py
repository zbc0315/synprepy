#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Time    : 2021/9/3 10:40
# @Author  : zhangbc0315@outlook.com
# @File    : spectators_predictor.py
# @Software: PyCharm

import torch
from torch.nn import functional as F
from torch_geometric.data import Batch

from nn_utils.gnn_data_utils import GnnDataUtils
from data_utils.tid_utils import TidUtils
from mcts_lite.mcts_path import MctsPath
from rxn_spectators_selector.ss_model import SSModel
from mcts_lite.mcts_config import MctsConfig


class SpectatorsPredictor:

    device = MctsConfig.ml_device
    basic_tids = TidUtils.load_basic_tids(fp=MctsPath.SPEC_TIDS_FP, sep='\n')
    length = len(basic_tids)
    cat_model = torch.load(MctsPath.CATS_MODEL_FP, map_location=torch.device(device)).to(device)
    sol_model = torch.load(MctsPath.SOLS_MODEL_FP, map_location=torch.device(device)).to(device)
    cat_model.eval()
    sol_model.eval()

    @classmethod
    def _get_global_vector(cls, tid: int):
        res = [0] * (cls.length + 1)
        if tid in cls.basic_tids:
            res[cls.basic_tids.index(tid)] = 1
        else:
            res[-1] = 1
        return [res]

    @classmethod
    def _get_p_batch(cls, p_data, num, global_vector_list):
        p_data_list = []
        for i in range(num):
            p_data_with_g = p_data.clone()
            p_data_with_g.g = global_vector_list[i]
            p_data_list.append(p_data_with_g)
        # p_data_list = [p_data.clone() for i in range(num)]
        return Batch.from_data_list(p_data_list)

    @classmethod
    def _add_global_to_reactants(cls, r_data_list, global_vector_list):
        for n, r_data in enumerate(r_data_list):
            r_data.g = global_vector_list[n]
        return r_data_list

    @classmethod
    def _get_gnn_data(cls, rd_reactants_list, rd_product, tids):
        global_vector_list = [torch.tensor(cls._get_global_vector(tid), dtype=torch.float) for tid in tids]

        r_data_list = [GnnDataUtils.gnn_data_from_mols(rd_reactants, 'cat_sol') for rd_reactants in rd_reactants_list]
        r_data_list = cls._add_global_to_reactants(r_data_list, global_vector_list)
        r_data_batch = Batch.from_data_list(r_data_list)

        p_data = GnnDataUtils.gnn_data_from_mol(rd_product, 'cat_sol')
        p_data_batch = cls._get_p_batch(p_data, len(r_data_list), global_vector_list)
        return r_data_batch, p_data_batch

    @classmethod
    def _batch_predict(cls, r_data_batch, p_data_batch, model):
        pred = model(r_data_batch, p_data_batch)
        # pred = F.softmax(pred)
        topk_probs, topk_idxes = torch.topk(pred, k=8, dim=1)
        return topk_idxes.tolist()

    @classmethod
    def batch_predict_cats_sols(cls, rd_reactants_list, rd_product, tids):
        r_data_batch, p_data_batch = cls._get_gnn_data(rd_reactants_list, rd_product, tids)
        r_data_batch = r_data_batch.to(cls.device)
        p_data_batch = p_data_batch.to(cls.device)
        cats_list = cls._batch_predict(r_data_batch, p_data_batch, cls.cat_model)
        sols_list = cls._batch_predict(r_data_batch, p_data_batch, cls.sol_model)
        return cats_list, sols_list


if __name__ == "__main__":
    reactants = ['CC(=O)O', 'CCO']
    products = ['CC(=O)OCC']
