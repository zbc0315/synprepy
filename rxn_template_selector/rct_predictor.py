#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Time    : 2022/2/20 11:05
# @Author  : zhangbc0315@outlook.com
# @File    : rct_predictor.py
# @Software: PyCharm

import logging
logging.basicConfig(level=logging.INFO)


import torch

from rxn_template_selector.rct_ret_dataset import RCTRETDataset
from rxn_template_selector.gnn_model import GNNModel
from data_utils import RxnTemplateData
from config.config import RCTSConfig
from utils.model_utils import ModelUtils


class RCTPredictor:

    def __init__(self, rct_config: RCTSConfig):
        self._rct_config = rct_config
        self._rxn_template_data = RxnTemplateData(rct_config.rxn_template_tsv_fp)
        best_model_state = ModelUtils.load_best_model(rct_config.model_dp)
        self._model = GNNModel(10, 5, 529+766+4, 15829).to(rct_config.device)
        self._model.load_state_dict(best_model_state)
        self._model.eval()

    def predict(self, data):
        with torch.no_grad():
            g, pred = self._model(data)
        return g, pred

    def predict_to_rt_smiles(self, data, topk: int) -> [str]:
        data = data.to(self._rct_config.device)
        _, pred_probs = self.predict(data)
        _, pred_idxes = pred_probs.topk(topk, 1, True, True)
        rt_smiles_list = [self._rxn_template_data.get_rt_smiles_by_ohi(idx) for idx in pred_idxes[0]]
        return rt_smiles_list


if __name__ == "__main__":
    pass
