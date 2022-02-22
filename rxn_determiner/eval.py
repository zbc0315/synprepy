#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Time    : 2021/9/2 10:27
# @Author  : zhangbc0315@outlook.com
# @File    : eval.py
# @Software: PyCharm

import torch
from tqdm import tqdm
from sklearn.metrics import roc_curve, auc
import matplotlib.pyplot as plt

from rxn_determiner.gcn_model import GcnModel
from rxn_determiner.dataloader_gcn import DataloaderGcn


torch.cuda.empty_cache()


class Eval:

    device = 'cuda'

    test_dataloader = DataloaderGcn('test')
    model = GcnModel(6, 12, 650).to(device)
    # model = GcnModel(6, 12, 0).to(device)
    model_state, _ = GcnModel.get_model_and_opt_state(True)
    model.load_state_dict(model_state)

    @classmethod
    def _get_best_threshold(cls, fpr, tpr, thresholds):
        ks_max = tpr[0] - fpr[0]
        best_thr = thresholds[0]
        idx = 0
        for i, thr, in enumerate(thresholds):
            if tpr[i] - fpr[i] > ks_max:
                ks_max = tpr[i] - fpr[i]
                best_thr = thresholds[i]
                idx = i
        return best_thr, idx

    @classmethod
    def draw_aoc(cls, pred, real):
        fpr, tpr, threshold = roc_curve(real, pred)
        roc_auc = auc(fpr, tpr)
        best_thr, i = cls._get_best_threshold(fpr, tpr, threshold)

        plt.figure()
        plt.figure(figsize=(5, 5))
        plt.plot(fpr, tpr, color='navy', lw=2, label=f'ROC curve (area={roc_auc})')
        plt.plot([0, 1], [0, 1], color='darkorange', lw=2, linestyle='--')
        plt.plot([fpr[i]], tpr[i], color='navy', marker='o', markerfacecolor='w')
        plt.xlim([0.0, 1.0])
        plt.ylim([0.0, 1.05])
        plt.xlabel('False Positive Rate')
        plt.ylabel('True Positive Rate')
        plt.show()
        print(f'auc: {roc_auc}')
        print(f'best thr: {best_thr}, fpr: {fpr[i]}, tpr: {tpr[i]}')

    @classmethod
    def eval_aoc(cls):
        total_pred = None
        total_real = None
        with tqdm(total=len(cls.test_dataloader))as pbar:
            for n, (r_data, p_data) in enumerate(cls.test_dataloader):
                r_data = r_data.to(cls.device)
                p_data = p_data.to(cls.device)

                pred = cls.model(r_data, p_data)
                real = r_data.y
                if total_pred is None:
                    total_pred = pred.to('cpu').tolist()
                    total_real = real.to('cpu').tolist()
                else:
                    total_pred.extend(pred.to('cpu').tolist())
                    total_real.extend(real.to('cpu').tolist())
                    # total_pred = torch.cat([total_pred, pred.to('cpu')])
                    # total_real = torch.cat([total_real, real.to('cpu')])
                pbar.update(1)

                if n >= len(cls.test_dataloader) - 1:
                    break
        cls.draw_aoc(total_pred, total_real)


if __name__ == "__main__":
    Eval.eval_aoc()
