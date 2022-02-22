#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Time    : 2022/1/13 21:41
# @Author  : zhangbc0315@outlook.com
# @File    : topk_acc_calculator.py
# @Software: PyCharm

from typing import Optional

import torch


class TopkAccCalculator:

    def __init__(self, max_k):
        self._max_k = max_k
        self._pred_idxes_list: Optional[torch.Tensor] = None
        self._pred_probs_list: Optional[torch.Tensor] = None
        self._real_idxes: Optional[torch.Tensor] = None

    def extend_pred_and_real(self, pred_probs_list, real_idxes):
        # pred_probs_list = pred_probs_list.to("cpu")
        # real_idxes = real_idxes.to("cpu")
        _, pred_idxes_list = pred_probs_list.topk(self._max_k, 1, True, True)
        if self._pred_idxes_list is None:
            self._pred_idxes_list = pred_idxes_list
            self._real_idxes = real_idxes
        else:
            self._pred_idxes_list = torch.cat((self._pred_idxes_list, pred_idxes_list))
            self._real_idxes = torch.cat((self._real_idxes, real_idxes))

    def get_topk_acc(self, k: int) -> float:
        num_tot = self._real_idxes.size()[0]
        # _, topk_idxes_list = self._pred_probs_list.topk(k, 1, True, True)
        topk_idxes_list = self._pred_idxes_list.t()
        correct = topk_idxes_list.eq(self._real_idxes.reshape((1, num_tot)).expand_as(topk_idxes_list))
        correct_k = correct[:k].reshape(-1).float().sum(0, keepdim=True)
        return correct_k[0].item() / num_tot


if __name__ == "__main__":
    pass
