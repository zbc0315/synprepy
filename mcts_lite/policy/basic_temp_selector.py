#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Time    : 2021/9/1 16:24
# @Author  : zhangbc0315@outlook.com
# @File    : basic_temp_selector.py
# @Software: PyCharm

import torch
import torch.nn.functional as F
from torch_geometric.data import Batch
from rdkit.Chem import AllChem

from mcts_lite.mcts_path import MctsPath
from rxn_template_selector.rts_gcn_model import RtsGcnModel
from nn_utils.gnn_data_utils import GnnDataUtils
from data_utils.tid_utils import TidUtils
from mcts_lite.mcts_config import MctsConfig


class BasicTempSelector:

    device = MctsConfig.ml_device
    model = RtsGcnModel(6, 12, 21715).to(device)
    model.load_state_dict(torch.load(MctsPath.MODEL_BT_FP, map_location=torch.device(device)))
    model.eval()
    basic_tids = TidUtils.load_basic_tids(fp=MctsPath.BASIC_TIDS_FP)

    @classmethod
    def predict(cls, rdmol):
        data = GnnDataUtils.gnn_data_from_mol(rdmol).to(cls.device)
        _, pred_batch = cls.model(Batch.from_data_list([data]))
        pred_batch = torch.softmax(pred_batch, dim=1)

        tid_probs = [(cls.basic_tids[i], prob.item()) for i, prob in enumerate(pred_batch[0])]
        sorted_tid_probs = sorted(tid_probs, key=lambda x: x[1], reverse=True)
        return sorted_tid_probs[:20]


if __name__ == "__main__":
    m = AllChem.MolFromSmiles('CCC')
    res = BasicTempSelector.predict(m)
    print(res[:10])
