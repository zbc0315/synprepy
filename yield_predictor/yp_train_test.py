#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Time    : 2022/2/12 22:05
# @Author  : zhangbc0315@outlook.com
# @File    : yp_train_test.py
# @Software: PyCharm

import logging

from tqdm import tqdm
import torch
from torch.nn import functional as F
from torch.optim import lr_scheduler
from torch_geometric.data import DataLoader
from matplotlib import pyplot as plt

from config import Config, YPConfig
from utils.model_utils import ModelUtils
from yield_predictor.yp_dataset import YPDataset
# from yield_predictor.yp_gnn_model import YPGNNModel as YPModel
from yield_predictor.yp_nn_model import YPNNModel as YPModel
from rxn_template_selector.rct_predictor import RCTPredictor

logging.basicConfig(level=logging.INFO)


class YPTrainTest:

    def __init__(self, config: Config):
        yp_config = config.yp_config
        batch_size = yp_config.batch_size
        self._epoch_num = yp_config.epoch_num
        lr_start = yp_config.lr_start
        lr_end = yp_config.lr_end
        self._device = yp_config.device
        self._model_dp = yp_config.model_dp

        # keys = ["px", "p_edge_attr", "p_edge_index"]
        _p_train_dataset = YPDataset("train", yp_config, 'p')
        self._p_train_dataloader = DataLoader(_p_train_dataset, batch_size=batch_size)
        _r_train_dataset = YPDataset("train", yp_config, 'r')
        self._r_train_dataloader = DataLoader(_r_train_dataset, batch_size=batch_size)

        _p_test_dataset = YPDataset("test", yp_config, 'p')
        self._p_test_dataloader = DataLoader(_p_test_dataset, batch_size=batch_size)
        _r_test_dataset = YPDataset("test", yp_config, 'r')
        self._r_test_dataloader = DataLoader(_r_test_dataset, batch_size=batch_size)

        self._model = YPModel(256, _r_train_dataset.num_global_features).to(self._device)
        # self._model = YPModel(_r_train_dataset.num_node_features, _r_train_dataset.num_edge_features, _r_train_dataset.num_global_features).to(self._device)
        self._model.load_state_dict(ModelUtils.load_best_model(config.yp_config.model_dp))
        self._optimizer = torch.optim.Adadelta(self._model.parameters(), lr=lr_start)
        self._scheduler = lr_scheduler.StepLR(self._optimizer, step_size=10, gamma=0.5, )
        self._reactant_vector_predictor = RCTPredictor(config.rcts_config)
        self._product_vector_predictor = RCTPredictor(config.retro_rcts_config)

    def train(self):
        for epoch in range(self._epoch_num):
            train_loss = self._train_epoch(epoch)
            test_loss = self.test(epoch)
            self._scheduler.step()
            logging.info(f"\ntrain loss: {train_loss}, test loss: {test_loss}")
            ModelUtils.save_model(self._model_dp, self._model, epoch, train_loss, 0, test_loss, 0)

    def _train_epoch(self, epoch: int):
        self._model.train()
        tot_loss = 0

        with tqdm(total=len(self._r_train_dataloader))as pbar:
            pbar.set_description_str(f"train epoch: {epoch}")
            for n, (r_batch_data, p_batch_data) in enumerate(zip(self._r_train_dataloader, self._p_train_dataloader)):
                if r_batch_data.key != p_batch_data.key:
                    raise ValueError("reactants and products are unmatched")
                r_batch_data = r_batch_data.to(self._device)
                p_batch_data = p_batch_data.to(self._device)
                rv_batch_data, _ = self._reactant_vector_predictor.predict(r_batch_data)
                pv_batch_data, _ = self._product_vector_predictor.predict(p_batch_data)
                rv_batch_data = torch.cat([rv_batch_data, r_batch_data.g], dim=1)
                #
                pred = self._model(rv_batch_data, pv_batch_data)
                # pred = self._model(r_batch_data, p_batch_data)
                real = p_batch_data.y
                loss = F.smooth_l1_loss(pred, real.reshape([real.size()[0], 1]))
                self._optimizer.zero_grad()
                loss.backward()
                self._optimizer.step()
                # for name, params in self._model.named_parameters():
                #     print('-->name: ', name, '-->grad_requires: ', params.requires_grad, '-->grad_value: ', params.grad)

                tot_loss += loss.item()
                pbar.set_postfix_str(f"loss: {tot_loss / (n + 1)}")
                pbar.update(1)
        return tot_loss / len(self._r_train_dataloader)

    def test(self, epoch: int):
        self._model.eval()
        tot_loss = 0
        pred_all = torch.tensor([], dtype=torch.float)
        real_all = torch.tensor([], dtype=torch.float)
        with torch.no_grad():
            with tqdm(total=len(self._r_test_dataloader))as pbar:
                pbar.set_description(f"test epoch: {epoch}")
                for n, (r_batch_data, p_batch_data) in enumerate(zip(self._r_test_dataloader, self._p_test_dataloader)):
                    if r_batch_data.key != p_batch_data.key:
                        raise ValueError("reactants and products are unmatched")
                    r_batch_data = r_batch_data.to(self._device)
                    p_batch_data = p_batch_data.to(self._device)

                    rv_batch_data, _ = self._reactant_vector_predictor.predict(r_batch_data)
                    pv_batch_data, _ = self._product_vector_predictor.predict(p_batch_data)
                    rv_batch_data = torch.cat([rv_batch_data, r_batch_data.g], dim=1)
                    #
                    pred = self._model(rv_batch_data, pv_batch_data)
                    # pred = self._model(r_batch_data, p_batch_data)
                    real = p_batch_data.y
                    loss = F.smooth_l1_loss(pred, real.reshape([real.size()[0], 1]))
                    tot_loss += loss.item()
                    pbar.set_postfix_str(f"loss: {tot_loss / (n + 1)}")
                    pbar.update(1)
                    pred_all = torch.cat([pred.to('cpu'), pred_all])
                    real_all = torch.cat([real.to('cpu'), real_all])
        plt.scatter(real_all.tolist(), pred_all.tolist(), alpha=0.3, s=2)
        plt.plot([0, 100], [0, 100])
        plt.show()
        return tot_loss / len(self._r_test_dataloader)


if __name__ == "__main__":
    proj_config = Config("../config.json")
    YPTrainTest(proj_config).train()
