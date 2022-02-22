#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Time    : 2021/8/31 12:59
# @Author  : zhangbc0315@outlook.com
# @File    : train_test.py
# @Software: PyCharm
import os

import torch
import torch.nn as nn
from tqdm import tqdm

from rxn_determiner.dataloader_gcn import DataloaderGcn
from rxn_determiner.gcn_model import GcnModel
from nn_utils.acc_utils import AccUtils
from nn_utils.param_utils import ParamUtils
from config import RDPath


class TrainTest:

    _EPOCH_START = 0
    _NUM_EPOCH = 200
    _NUM_BATCH = 1000
    _DEVICE = 'cuda'
    _LR_START = 0.0005
    _LR_END = 1e-6

    _train_data_loader = DataloaderGcn('train')
    _test_data_loader = DataloaderGcn('test')
    _data_param = DataloaderGcn.get_data_param()
    _num_global = 0
    # _num_global = _data_param['num global features']
    _model = GcnModel(num_node_features=_data_param['num node features'],
                      num_edge_features=_data_param['num edge features'],
                      num_global_features=_num_global).to(_DEVICE)
    # _model = GcnModel(num_node_features=_data_param['num node features'],
    #                   num_edge_features=_data_param['num edge features'],
    #                   num_global_features=0).to(_DEVICE)
    _optimizer = torch.optim.Adagrad(_model.parameters(), lr=_LR_START)
    _loss_func = nn.BCELoss()

    # _model_opt_fp = "C:\\Users\\zhang\\Documents\\Data\\RxnMolPredictor\\Mcts\\Models\\rd_model_opt_122_a0.9700794741305372_p0.9702879087815472_r0.9631502495135691_.pt"
    # _model_opt_state = torch.load(_model_opt_fp)
    # _model.load_state_dict(_model_opt_state['model'])

    @classmethod
    def train(cls):
        cls._model.train()
        total_loss = 0
        total_tp = 0
        total_tn = 0
        total_fp = 0
        total_fn = 0
        with tqdm(total=len(cls._train_data_loader))as pbar:
            for n, (r_data, p_data) in enumerate(cls._train_data_loader):
                r_data = r_data.to(cls._DEVICE)
                p_data = p_data.to(cls._DEVICE)
                pred = cls._model(r_data, p_data)
                real = torch.tensor(r_data.y, dtype=torch.float)
                real = real.to(cls._DEVICE)
                loss = cls._loss_func(pred, real)
                loss.backward()
                cls._optimizer.step()

                total_loss += loss.item()
                tp, tn, fp, fn = AccUtils.calc_tfpn(pred, real)

                total_tp += tp
                total_tn += tn
                total_fp += fp
                total_fn += fn
                acc, pre, rec = AccUtils.calc_acc_pre_rec(total_tp, total_tn, total_fp, total_fn)

                pbar.set_postfix_str(f'loss: {total_loss / (n + 1)}, acc: {acc}, pre: {pre}, recall: {rec}, '
                                     f'tp: {total_tp}, tn: {total_tn}, fp: {total_fp}, fn: {total_fn}')
                pbar.update(1)
                if n >= len(cls._train_data_loader) - 1:
                    break
        acc, pre, rec = AccUtils.calc_acc_pre_rec(total_tp, total_tn, total_fp, total_fn)
        return acc, pre, rec

    @classmethod
    def test(cls):
        cls._model.eval()
        total_tp = 0
        total_tn = 0
        total_fp = 0
        total_fn = 0
        with tqdm(total=len(cls._test_data_loader))as pbar:
            for n, (r_data, p_data) in enumerate(cls._test_data_loader):
                if n in [50]:
                    continue
                r_data = r_data.to(cls._DEVICE)
                p_data = p_data.to(cls._DEVICE)
                pred = cls._model(r_data, p_data)
                real = torch.tensor(r_data.y, dtype=torch.float).to(cls._DEVICE)

                tp, tn, fp, fn = AccUtils.calc_tfpn(pred, real)
                total_tp += tp
                total_tn += tn
                total_fp += fp
                total_fn += fn
                acc, pre, rec = AccUtils.calc_acc_pre_rec(total_tp, total_tn, total_fp, total_fn)

                pbar.set_postfix_str(f'acc: {acc}, pre: {pre}, recall: {rec}, '
                                     f'tp: {total_tp}, tn: {total_tn}, fp: {total_fp}, fn: {total_fn}')
                pbar.update(1)
                if n >= len(cls._test_data_loader) - 1:
                    break
        acc, pre, rec = AccUtils.calc_acc_pre_rec(total_tp, total_tn, total_fp, total_fn)
        print(f'acc: {acc}, pre: {pre}, recall: {rec}, '
              f'tp: {total_tp}, tn: {total_tn}, fp: {total_fp}, fn: {total_fn}')

    @classmethod
    def _save_model(cls, e, acc, pre, rec, dp):
        model_fn = f'model_opt_{e}_a{acc}_p{pre}_r{rec}_.pt'
        torch.save({'model': cls._model.state_dict(),
                    'opt': cls._optimizer.state_dict()},
                   os.path.join(dp, model_fn))

    @classmethod
    def _load_model(cls, dp):
        fns = sorted(os.listdir(dp), key=lambda x: float(x.split('_')[3][1:]))
        best_fn = fns[-1]
        best_fp = os.path.join(dp, best_fn)
        current_e = int(best_fn.split('_')[2])
        print(f'load model: {best_fp}')
        state = torch.load(best_fp)
        cls._model.load_state_dict(state['model'])
        cls._optimizer.load_state_dict(state['opt'])
        return current_e

    @classmethod
    def train_and_test(cls):
        dp = RDPath.MODEL_CAT_SOL_DP if cls._num_global != 0 else RDPath.MODEL_NO_CAT_SOL_DP
        current_e = 0
        current_e = cls._load_model(dp)
        # ParamUtils.adjust_lr(cls._LR_START, cls._LR_END, cls._optimizer, current_e)
        for e in range(max(cls._EPOCH_START, current_e), cls._NUM_EPOCH):

            acc, pre, rec = cls.train()
            # cls.test()
            cls._save_model(e, acc, pre, rec, dp)
            ParamUtils.adjust_lr(cls._LR_START, cls._LR_END, cls._optimizer, e)


if __name__ == "__main__":
    TrainTest.train_and_test()
    # TrainTest.test()
