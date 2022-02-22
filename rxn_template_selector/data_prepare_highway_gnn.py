#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Time    : 2021/8/16 15:05
# @Author  : zhangbc0315@outlook.com
# @File    : data_prepare_highway_gnn.py
# @Software: PyCharm

import os

import torch

from config import BTPath
from rxn_template_selector.rts_gcn_model import RtsGcnModel


class DataPrepareHighwayGnn:

    @classmethod
    def _get_gcn_model(cls, gcn_models_dp: str):
        # gcn_model = RtsGcnModel(6, 12, 21715)
        gcn_model = RtsGcnModel(6, 12, 2133)
        gcn_model_fp = RtsGcnModel.get_model_fp(gcn_models_dp, "best")
        gcn_model_fp = 'C:\\Users\\zhang\\Documents\\Data\\RxnMolPredictor\\Mcts\\Models\\' \
                       'bt_model_200_2.185096339638548_0.39687987110802564.pt'
        model_state = torch.load(gcn_model_fp)
        gcn_model.load_state_dict(model_state)
        gcn_model.eval()
        return gcn_model

    @classmethod
    def prepare_data_fp(cls, source_fp: str, target_fp: str, gcn_model):
        data = torch.load(source_fp)
        g, _ = gcn_model(data)
        out = {"x": g, "y": data.y}
        torch.save(out, target_fp)

    @classmethod
    def prepare_data_dp(cls, source_dp: str, target_dp: str, gcn_model):
        for n, source_fn in enumerate(os.listdir(source_dp)):
            if os.path.exists(os.path.join(target_dp, source_fn)):
                print(f"{n} - parsed: {source_fn}")
            else:
                print(f"{n} - parsing: {source_fn}")
                cls.prepare_data_fp(os.path.join(source_dp, source_fn),
                                    os.path.join(target_dp, source_fn),
                                    gcn_model)

    @classmethod
    def process(cls, gcn_models_dp: str, source_train_dp: str, source_test_dp: str,
                target_train_dp: str, target_test_dp: str):
        gcn_model = cls._get_gcn_model(gcn_models_dp)
        cls.prepare_data_dp(source_train_dp, target_train_dp, gcn_model)
        cls.prepare_data_dp(source_test_dp, target_test_dp, gcn_model)


if __name__ == "__main__":
    source_data_dp = "C:\\Users\\zhang\\Documents\\Data\\RxnPredictor\\GNN\\BaseTemplate_AtomicNum"
    # DataPrepareHighwayGnn.process(os.path.join(source_data_dp, "model_neg"),
    #                               os.path.join(source_data_dp, "train_processed_batch_1000"),
    #                               os.path.join(source_data_dp, "test_processed_batch_1000"),
    #                               BTPath.TRAIN_BATCH_DATA_DP,
    #                               BTPath.TEST_BATCH_DATA_DP)

    source_comp_data_dp = "C:\\Users\\zhang\\Documents\\Data\\RxnMolPredictor\\CompTemp\\"
    DataPrepareHighwayGnn.process(os.path.join(source_data_dp, "model_neg"),
                                  os.path.join(source_comp_data_dp, 'gcn\\train_batch_1000_data'),
                                  os.path.join(source_comp_data_dp, 'gcn\\test_batch_1000_data'),
                                  os.path.join(source_comp_data_dp, 'highway_gcn\\train_batch_1000_data'),
                                  os.path.join(source_comp_data_dp, 'highway_gcn\\test_batch_1000_data'))
