#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Time    : 2021/8/24 14:30
# @Author  : zhangbc0315@outlook.com
# @File    : train_test_gcn.py
# @Software: PyCharm
import os

import torch.nn.functional as F
import torch.optim
from tqdm import tqdm
import pandas as pd

from nn_utils import AccUtils
from nn_utils import ParamUtils
from nn_utils import ModelNameUtils
from rxn_template_selector.dataloader_gcn import DataloaderGcn
from rxn_template_selector.rts_gcn_model import RtsGcnModel
from config import BTPath, CTPath
from mcts_lite.mcts_path import MctsPath


class TrainTestGcn:

    _EPOCH_START = 0
    _NUM_EPOCH = 200
    _NUM_BATCH = 1000
    _DEVICE = "cuda"
    _LR_START = 0.0005
    _LR_END = 1e-7

    _proj = 'base_temp'
    _train_data_loader = DataloaderGcn(_proj, 'train')
    _test_data_loader = DataloaderGcn(_proj, 'test')
    _data_param = _test_data_loader.get_data_param()

    _model = RtsGcnModel(num_node_features=_data_param['num node features'],
                         num_edge_features=_data_param['num edge features'],
                         num_out_channels=_data_param['num out channels']).to(_DEVICE)
    _optimizer = torch.optim.Adagrad(_model.parameters(), lr=_LR_START)
    # _optimizer = torch.optim.Adam(_model.parameters(), lr=_LR_START)

    _model_dps = {'base_temp': BTPath.CBE_MODEL_DP,
                  'comp_temp': CTPath.MODEL_DP}
    # _best_model_fp = MctsPath.MODEL_BT_FP
    # model_state = torch.load(_best_model_fp)
    # _model.load_state_dict(model_state)
    _best_model_fp = RtsGcnModel.get_model_fp(_model_dps[_proj], sort_key='latest')
    model_opt = torch.load(_best_model_fp)
    _model.load_state_dict(model_opt)

    @classmethod
    def _init_model_and_opt(cls, init_epoch: bool):
        model_state, opt_state, epoch = RtsGcnModel.get_model_and_opt(cls._model_dps[cls._proj])
        if init_epoch:
            cls._EPOCH_START = epoch
        if model_state is not None:
            cls._model.load_state_dict(model_state)
            cls._optimizer.load_state_dict(opt_state)

    @classmethod
    def parse_y(cls, ys):
        return ys.max(dim=1)[1]

    @classmethod
    def add_topk_acc(cls, total_topk_acc, topk_acc):
        if len(total_topk_acc["num"]) == 0:
            total_topk_acc["num"] = topk_acc["num"]
            total_topk_acc["acc"] = []
            for i in topk_acc["num"]:
                total_topk_acc["acc"].append([])
        for i, acc in enumerate(topk_acc["acc"]):
            total_topk_acc["acc"][i].append(acc)
        return total_topk_acc

    @classmethod
    def calc_mean_topk(cls, total_topk_acc, data_nums):
        res = {"num": total_topk_acc["num"],
               "acc": [0]*len(total_topk_acc["num"])}
        sum_data_num = sum(data_nums)
        for i, accs in enumerate(total_topk_acc["acc"]):
            for acc, data_num in zip(accs, data_nums):
                res["acc"][i] += acc * data_num
            res["acc"][i] = res["acc"][i] / sum_data_num
        return res

    @classmethod
    def train(cls):
        cls._model.train()
        total_loss = 0
        total_correct = 0
        with tqdm(total=len(cls._train_data_loader))as pbar:
            for n, data in enumerate(cls._train_data_loader):
                data = data.to(cls._DEVICE)
                _, pred = cls._model(data)
                y = cls._model.change_y(data.y)
                loss = F.cross_entropy(pred, cls.parse_y(y))
                loss.backward()
                cls._optimizer.step()

                total_loss += loss.item()
                total_correct += cls.parse_y(pred).eq(cls.parse_y(y)).sum().item()
                pbar.set_postfix_str(f"loss: {total_loss / (n + 1)}; acc: {total_correct / ((n + 1) * cls._NUM_BATCH)}")
                pbar.update(1)
        return total_loss / len(cls._train_data_loader), total_correct / (len(cls._train_data_loader) * cls._NUM_BATCH)

    @classmethod
    def _get_short_acc(cls, acc):
        res = {}
        for k in acc.keys():
            res[k] = [acc[k][0], acc[k][-1]]
        return res

    @classmethod
    def test(cls, ks: [int] = None):
        ks = ks if ks is not None else [1, 5, 10, 30, 50]
        cls._model.eval()
        total_topk_acc = {"num": [],
                          "acc": []}
        data_nums = []
        with tqdm(total=len(cls._test_data_loader))as pbar:
            for n, data in enumerate(cls._test_data_loader):
                data = data.to(cls._DEVICE)
                _, pred = cls._model(data)
                y = cls._model.change_y(data.y)
                topk_acc, data_num = AccUtils.get_topk_acc(pred, y, ks=ks)
                total_topk_acc = cls.add_topk_acc(total_topk_acc, topk_acc)
                data_nums.append(data_num)
                mean_topk_acc = cls.calc_mean_topk(total_topk_acc.copy(), data_nums)
                pbar.set_postfix_str(f"{cls._get_short_acc(mean_topk_acc)}")
                pbar.update(1)
        return mean_topk_acc

    @classmethod
    def save_model_and_opt(cls, epoch, loss, acc):
        model_fn = ModelNameUtils.get_model_name(epoch, loss, acc)
        model_fp = os.path.join(cls._model_dps[cls._proj], model_fn)
        res = {'model': cls._model.state_dict(), 'opt': cls._optimizer.state_dict()}
        torch.save(res, model_fp)

    @classmethod
    def train_test(cls):
        # cls._init_model_and_opt(init_epoch=True)
        for e in range(cls._EPOCH_START, cls._NUM_EPOCH):
            ParamUtils.adjust_lr(cls._LR_START, cls._LR_END, cls._optimizer, e)
            loss, acc = cls.train()
            cls.test()
            cls.save_model_and_opt(e, loss, acc)

    @classmethod
    def only_test(cls):
        mean_topk_acc = cls.test(list(range(1, 51)))
        df = pd.DataFrame(mean_topk_acc)
        df.to_csv(CTPath.GCN_RESULT_FP, sep='\t', encoding='utf-8')


if __name__ == "__main__":
    TrainTestGcn.only_test()
    # TrainTestGcn.train_test()
