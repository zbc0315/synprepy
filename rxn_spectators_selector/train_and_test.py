# -*- coding: utf-8 -*-
# @Time     : 2021/4/28 16:17
# @Author   : zbc
# @Email    : zhangbc0315@outlook.com
# @FileName : train_and_test.py
# @Software : PyCharm

import os

import torch
import torch.nn.functional as F
from torch.utils.tensorboard import SummaryWriter
from tqdm import tqdm

from rxn_spectators_selector.rxn_data import RxnData
from rxn_spectators_selector.ss_model import SSModel
# from rxn_spectators_selector.ss_pool_gcn_model import SSPoolGcnModel as SSModel
from config import SSPath


class TrainAndTest:

    _num_epochs = 300
    _device = "cuda"
    _batch = 500
    _model_dps = {"cat": SSPath.CAT_GNN_MODEL_DP,
                  "sol": SSPath.SOL_GNN_MODEL_DP}

    @classmethod
    def _parse_y(cls, y):
        return y.max(dim=1)[1]

    @classmethod
    def train(cls, model, train_set, optimizer):
        model.train()
        loss_all = 0
        correct_all = 0
        optimizer.zero_grad()
        with tqdm(total=int(train_set.len()/cls._batch)) as pbar:
            for n, (r_batch_data, p_batch_data) in enumerate(train_set.get_batch_data()):
                r_batch_data = r_batch_data.to(cls._device)
                p_batch_data = p_batch_data.to(cls._device)
                output = model(r_batch_data, p_batch_data)
                loss = F.nll_loss(output, cls._parse_y(r_batch_data.y))
                # loss = F.nll_loss(output, r_batch_data.y)
                # loss = F.nll_loss(output, r_batch_data.y)
                loss.backward()
                loss_all += r_batch_data.num_graphs * loss.item()
                optimizer.step()
                pred = cls._parse_y(output)
                correct = pred.eq(cls._parse_y(r_batch_data.y)).sum().item()
                correct_all += correct
                pbar.set_postfix_str(f"loss: {loss_all / ((n+1)*cls._batch)} "
                                     f"acc: {correct_all / ((n+1)*cls._batch)}")
                pbar.update(1)
        return loss_all / train_set.len()

    @classmethod
    def test(cls, model, test_set):
        model.eval()
        correct = 0
        loss = 0
        with tqdm(total=int(test_set.len() / cls._batch)) as pbar:
            for n, (r_batch_data, p_batch_data) in enumerate(test_set.get_batch_data()):
                r_batch_data = r_batch_data.to(cls._device)
                p_batch_data = p_batch_data.to(cls._device)
                out_put = model(r_batch_data, p_batch_data)
                pred = cls._parse_y(out_put)
                correct += pred.eq(cls._parse_y(r_batch_data.y)).sum().item()
                loss += F.nll_loss(out_put, cls._parse_y(r_batch_data.y)).item()
                pbar.update(1)
        return loss / test_set.len(), correct / test_set.len()

    @classmethod
    def _get_model_fp(cls, model_dp, e, test_loss, test_correct):
        return os.path.join(model_dp, f"model_{e}_l{test_loss}_c{test_correct}.pt")

    @classmethod
    def train_and_test(cls, role):
        writer = SummaryWriter(SSPath.TENSOR_BOARD_DP)
        model_dp = cls._model_dps[role]
        train_set = RxnData(role, "train")
        test_set = RxnData(role, "test")
        print(f"num node features: {train_set.num_node_features}\n"
              f"num edge features: {train_set.num_edge_features}\n"
              f"num g: {train_set.num_g}\n"
              f"out channels: {train_set.out_channels}")
        model = SSModel(train_set.num_node_features,
                        train_set.num_edge_features,
                        train_set.num_g,
                        train_set.out_channels).to(cls._device)
        # optimizer = torch.optim.Adam(model.parameters(), lr=0.0002)
        optimizer = torch.optim.Adagrad(model.parameters(), lr=0.0005)
        for e in range(cls._num_epochs):
            train_loss = cls.train(model, train_set, optimizer)
            writer.add_scalar(f"loss/train_{role}", train_loss, e)
            test_loss, test_correct = cls.test(model, test_set)
            writer.add_scalar(f"loss/test_{role}", test_loss, e)
            torch.save(model, cls._get_model_fp(model_dp, e, test_loss, test_correct))


if __name__ == '__main__':
    TrainAndTest.train_and_test("sol")
