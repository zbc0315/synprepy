#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Time    : 2022/1/11 10:30
# @Author  : zhangbc0315@outlook.com
# @File    : rct_ret_train_test.py
# @Software: PyCharm
import logging

from tqdm import tqdm
import torch
from torch.nn import functional as F
from torch.optim import lr_scheduler
from torch_geometric.data import DataLoader, Batch

from utils.torch_utils import TorchUtils
from config.config import Config
from utils.model_utils import ModelUtils
from rxn_template_selector.rct_ret_dataset import RCTRETDataset
from rxn_template_selector.gnn_model import GNNModel
from data_utils.rxn_template_data import RxnTemplateData
from utils.topk_acc_calculator import TopkAccCalculator

logging.basicConfig(level=logging.INFO)


class RCTRETTrainTest:

    def __init__(self, config):
        batch_size = config.batch_size
        self._epoch_num = config.epoch_num
        lr_start = config.lr_start
        lr_end = config.lr_end
        self._device = config.device
        self._model_dp = config.model_dp

        _train_dataset = RCTRETDataset("train", config)
        self._train_dataloader = DataLoader(_train_dataset, batch_size=batch_size)

        _test_dataset = RCTRETDataset("test", config)
        self._test_dataloader = DataLoader(_test_dataset, batch_size=batch_size)

        self._rxn_template_data = RxnTemplateData(config.rxn_template_tsv_fp)
        self._rxn_template_data.init_filter(config.min_num_covered_rxns_by_rxn_template)

        self._model = GNNModel(_train_dataset.num_node_features,
                               _train_dataset.num_edge_features,
                               _train_dataset.num_global_features,
                               len(self._rxn_template_data)).to(self._device)
        self._optimizer = torch.optim.Adam(self._model.parameters(), lr=lr_start)
        self._scheduler = lr_scheduler.StepLR(self._optimizer, step_size=10, gamma=0.5, )
        # self._scheduler = lr_scheduler.ReduceLROnPlateau(self._optimizer, 'min', factor=0.5, min_lr=lr_end)

    def train(self):
        for epoch in range(self._epoch_num):
            train_loss, train_acc = self._train_epoch(epoch)
            test_loss, test_acc = self.test(epoch)
            self._scheduler.step()
            logging.info(f"\ntrain loss: {train_loss}, acc: {train_acc}\n"
                         f"test loss: {test_loss}, acc: {test_acc}")
            ModelUtils.save_model(self._model_dp, self._model, epoch,
                                  train_loss, train_acc,
                                  test_loss, test_acc)

    def _train_epoch(self, epoch: int):
        self._model.train()
        tot_loss = 0
        tot_correct = 0
        tot_mol = 0
        with tqdm(total=len(self._train_dataloader))as pbar:
            pbar.set_description_str(f"train epoch: {epoch}")
            for n, batch_data in enumerate(self._train_dataloader):
                self._optimizer.zero_grad()
                batch_data = batch_data.to(self._device)
                g, pred = self._model(batch_data)
                loss = F.cross_entropy(pred, batch_data.y)

                loss.backward()
                self._optimizer.step()

                num_mols = batch_data.num_graphs
                correct = pred.max(dim=1)[1].eq(batch_data.y).sum().item()
                tot_correct += correct
                tot_mol += num_mols
                tot_loss += (loss.item()*num_mols)
                pbar.set_postfix_str(f"loss: {tot_loss / tot_mol}, acc: {tot_correct / tot_mol}")
                pbar.update(1)
        return tot_loss / tot_mol, tot_correct / tot_mol

    def test(self, epoch: int):
        self._model.eval()
        tot_loss = 0
        tot_correct = 0
        tot_mol = 0
        tac = TopkAccCalculator(51)
        with torch.no_grad():
            with tqdm(total=len(self._test_dataloader))as pbar:
                pbar.set_description_str(f"test epoch: {epoch}")
                for n, batch_data in enumerate(self._test_dataloader):
                    batch_data = batch_data.to(self._device)
                    g, pred = self._model(batch_data)
                    loss = F.cross_entropy(pred, batch_data.y)

                    num_mols = batch_data.num_graphs
                    correct = pred.max(dim=1)[1].eq(batch_data.y).sum().item()
                    tot_correct += correct
                    tot_mol += num_mols
                    tot_loss += (loss.item() * num_mols)
                    pbar.set_postfix_str(f"loss: {tot_loss / tot_mol}, acc: {tot_correct / tot_mol}")
                    pbar.update(1)
                    tac.extend_pred_and_real(pred, batch_data.y)
        logging.info(f"test topk acc(1, 5, 10, 20, 50): "
                     f"{tac.get_topk_acc(1)}, "
                     f"{tac.get_topk_acc(5)}, "
                     f"{tac.get_topk_acc(10)}, "
                     f"{tac.get_topk_acc(20)}, "
                     f"{tac.get_topk_acc(50)}")
        return tot_loss / tot_mol, tot_correct / tot_mol

    def eval(self):
        self._model.load_state_dict(ModelUtils.load_best_model(self._model_dp))
        tac = TopkAccCalculator(60)
        self._model.eval()
        with torch.no_grad():
            with tqdm(total=len(self._test_dataloader))as pbar:
                pbar.set_description_str(f"eval")
                for n, batch_data in enumerate(self._test_dataloader):
                    batch_data = batch_data.to(self._device)
                    g, pred = self._model(batch_data)
                    tac.extend_pred_and_real(pred, batch_data.y)
                    pbar.update(1)
        for i in [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]:
            logging.info(f"top{i}: {tac.get_topk_acc(i)}")
        # logging.info(f"\n"
        #              f"top1: {tac.get_topk_acc(1)}\n"
        #              f"top2: {tac.get_topk_acc(2)}\n"
        #              f"top3: {tac.get_topk_acc(3)}\n"
        #              f"top4: {tac.get_topk_acc(4)}\n"
        #              f"top5: {tac.get_topk_acc(5)}")


if __name__ == "__main__":
    proj_config = Config("../config.json")
    # RCTRETTrainTest(proj_config.retro_rcts_config).train()
    RCTRETTrainTest(proj_config.rets_config).train()
    # RCTRETTrainTest(proj_config.rcts_config).eval()
