#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Time    : 2021/8/17 15:24
# @Author  : zhangbc0315@outlook.com
# @File    : acc_utils.py
# @Software: PyCharm

import torch


class AccUtils:

    @classmethod
    def calc_acc_pre_rec(cls, tp, tn, fp, fn):
        acc = (tp+tn) / (tp+tn+fp+fn)
        pre = tp/(tp+fp) if tp + fp > 0 else -1
        rec = tp/(tp+fn) if tp + fn > 0 else -1
        return acc, pre, rec

    @classmethod
    def calc_tfpn(cls, pred_list, real_list):
        tp = 0
        tn = 0
        fp = 0
        fn = 0
        threshold = 0.5
        for pred, real in zip(pred_list, real_list):
            if pred >= threshold and real == 1:
                tp += 1
            elif pred < threshold and real == 0:
                tn += 1
            elif pred >= threshold and real == 0:
                fp += 1
            else:
                fn += 1
        return tp, tn, fp, fn

    # ========================

    @classmethod
    def get_acc_idx(cls, pred_probs: torch.Tensor, real_probs: torch.Tensor, max_idx: int = 50):
        pred_idxes = torch.topk(pred_probs, k=50)[1]
        real_idx = real_probs.max(dim=1)[1]
        for n, one_pred_idxes in enumerate(pred_idxes):
            one_real_idx = real_idx[n]
            try:
                yield list(one_pred_idxes).index(one_real_idx.item())
            except Exception as e:
                yield 100

    @classmethod
    def get_topk_acc(cls, pred_probs: torch.Tensor, real_probs: torch.Tensor, max_idx: int = 50):
        correct_idxes = torch.tensor(list(cls.get_acc_idx(pred_probs, real_probs, max_idx=max_idx)))
        res = {'num': [], 'acc': []}
        for num in [1, 5, 10, 30, 50]:
            count = correct_idxes.gt(num-1).sum().item()
            acc = 1-count/len(correct_idxes)
            res['num'].append(num)
            res['acc'].append(acc)
        return res, len(correct_idxes)


if __name__ == "__main__":
    pass
