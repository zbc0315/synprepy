#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Time    : 2021/9/3 16:24
# @Author  : zhangbc0315@outlook.com
# @File    : comp_temp_selector.py
# @Software: PyCharm
import torch
from torch_geometric.data import Batch
from rdkit.Chem import AllChem

from rxn_template_selector.rts_highway_dnn_model import HighwayDnnModel
from rxn_template_selector.rts_gcn_model import RtsGcnModel
from nn_utils.gnn_data_utils import GnnDataUtils
from data_utils.tid_utils import TidUtils
from mcts_lite.mcts_path import MctsPath
from mcts_lite.mcts_config import MctsConfig


class CompTempSelector:

    device = MctsConfig.ml_device
    gnn_model = RtsGcnModel(6, 12, 21715).to(device)
    gnn_model.load_state_dict(torch.load(MctsPath.MODEL_BT_FP, map_location=torch.device(device)))
    gnn_model.eval()

    dnn_model = HighwayDnnModel(256, 2133).to(device)
    state = torch.load(MctsPath.MODEL_CT_FP, map_location=torch.device(device))
    dnn_model.load_state_dict(state['model'])
    dnn_model.eval()

    comp_tids = TidUtils.load_comp_tid_to_idx(fp=MctsPath.COMP_TIDS_FP, res_type='list')

    @classmethod
    def _batch_predict(cls, rdmol_list):
        gnn_data_list = [GnnDataUtils.gnn_data_from_mol(rdmol) for rdmol in rdmol_list]
        gnn_data_batch = Batch.from_data_list(gnn_data_list)
        gnn_data_batch = gnn_data_batch.to(cls.device)
        dnn_input, _ = cls.gnn_model(gnn_data_batch)
        pred = cls.dnn_model(dnn_input)
        topk_probs, topk_idxes = torch.topk(pred, k=10, dim=1)
        for idxes in topk_idxes:
            res = [cls.comp_tids[idx] for idx in idxes]
            yield res

    @classmethod
    def batch_predict(cls, rdmol_list):
        for n in range(len(rdmol_list) // 1000 + 1):
            start = n * 1000
            end = start + 1000
            end = min(len(rdmol_list), end)
            for tids in cls._batch_predict(rdmol_list[start: end]):
                yield tids

    @classmethod
    def batch_predict_by_inchis(cls, inchi_list):
        rdmol_list = [AllChem.MolFromInchi(inchi) for inchi in inchi_list]
        for res in cls.batch_predict(rdmol_list):
            yield res


def test():

    mol1 = AllChem.MolFromSmiles('CCC')
    mol2 = AllChem.MolFromSmiles('CC')
    for res in CompTempSelector.batch_predict([mol1, mol2]):
        print(res)


if __name__ == "__main__":
    test()
