#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Time    : 2021/8/16 16:16
# @Author  : zhangbc0315@outlook.com
# @File    : train_test.py
# @Software: PyCharm

import os

from tqdm import tqdm
import torch.nn.functional as F
import torch.optim.optimizer
import pandas as pd

from config import BTPath, CTPath
from nn_utils.acc_utils import AccUtils
from nn_utils.param_utils import ParamUtils
from rxn_template_selector.dataloader_highway import DataloaderHighway
from rxn_template_selector.rts_highway_dnn_model import HighwayDnnModel


class TrainTest:

    _epoch_num = 200
    _NUM_INPUT = 256
    _num_batch = 1000
    # _out_channel = 21715
    _out_channel = 2133

    _proj = 'comp_temp'
    _model_dps = {'comp_temp': 'C:\\Users\\zhang\\Documents\\Data\\RxnMolPredictor\\CompTemp\\highway_gcn\\models',
                  'base_temp': BTPath.HG_MODELS_DP}
    _model_dp = _model_dps[_proj]

    _device = "cpu"
    _model, _epoch_begin, _opt_state = HighwayDnnModel.load_model(_NUM_INPUT, _out_channel, _device, _model_dp)
    _lr_start = 0.0005
    _lr_end = 1e-5
    _optimizer = torch.optim.Adagrad(_model.parameters(), lr=_lr_start)
    # _optimizer = torch.optim.Adam(_model.parameters(), lr=_lr_start)
    if _opt_state is not None:
        _optimizer.load_state_dict(_opt_state)

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
    def train(cls, dataloader):
        cls._model.train()
        total_loss = 0
        total_correct = 0
        with tqdm(total=len(dataloader)) as pbar:
            for n, data in enumerate(dataloader):
                # if n >= 10:
                #     break
                x = data['x'].to(cls._device)
                y = data['y'].to(cls._device)
                # parsed_y = cls.parse_y(y)
                parsed_y = cls.parse_y(y)
                pred = cls._model(x)
                loss = F.cross_entropy(pred, parsed_y)
                loss.backward()
                cls._optimizer.step()

                total_loss += loss.item()
                total_correct += cls.parse_y(pred).eq(parsed_y).sum().item()
                pbar.set_postfix_str(f"loss: {total_loss / (n+1)}; acc: {total_correct / ((n + 1) * cls._num_batch)}")
                pbar.update(1)
        return total_loss / len(dataloader), total_correct / (len(dataloader) * cls._num_batch)

    @classmethod
    def test(cls, dataloader, ks: [int] = None):
        ks = ks if ks is not None else [1, 5, 10, 30, 50]
        cls._model.eval()
        total_correct = 0
        total_topk_acc = {"num": [],
                          "acc": []}
        data_nums = []
        with tqdm(total=len(dataloader)) as pbar:
            for n, data in enumerate(dataloader):
                x = data['x'].to(cls._device)
                y = data['y'].to(cls._device)
                pred = cls._model(x)
                topk_acc, data_num = AccUtils.get_topk_acc(pred, y, ks=ks)
                total_topk_acc = cls.add_topk_acc(total_topk_acc, topk_acc)
                data_nums.append(data_num)
                mean_topk_acc = cls.calc_mean_topk(total_topk_acc.copy(), data_nums)
                pbar.set_postfix_str(f"{mean_topk_acc}")
                pbar.update(1)
        return mean_topk_acc

    @classmethod
    def save_model_and_opt(cls, epoch, loss, acc):
        model_fp = os.path.join(cls._model_dp, f"{epoch}_{loss}_{acc}_.pt")
        res = {'model': cls._model.state_dict(), 'opt': cls._optimizer.state_dict()}
        torch.save(res, model_fp)

    @classmethod
    def train_test(cls):
        dataloader_train = DataloaderHighway(cls._proj, "train")
        dataloader_test = DataloaderHighway(cls._proj, "test")
        for e in range(cls._epoch_begin+1, cls._epoch_num):
            ParamUtils.adjust_lr(cls._lr_start, cls._lr_end, cls._optimizer, e)
            loss, acc = cls.train(dataloader_train)
            cls.test(dataloader_test)
            cls.save_model_and_opt(e, loss, acc)

    @classmethod
    def only_test(cls):
        dataloader_test = DataloaderHighway(cls._proj, "test")
        mean_topk_acc = cls.test(dataloader_test, ks=list(range(1, 51)))
        df = pd.DataFrame(mean_topk_acc)
        df.to_csv(CTPath.HIGHWAY_RESULT_FP, sep='\t', encoding='utf-8', index=False)
        # print(mean_topk_acc)


if __name__ == "__main__":
    # TrainTest.only_test()
    TrainTest.train_test()
