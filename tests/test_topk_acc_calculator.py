#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Time    : 2022/1/13 23:37
# @Author  : zhangbc0315@outlook.com
# @File    : test_topk_acc_calculator.py
# @Software: PyCharm
from unittest import TestCase

import torch

from utils.topk_acc_calculator import TopkAccCalculator


class TestTopkAccCalculator(TestCase):

    def test_top1_acc(self):
        pred_probs_list = torch.tensor([[0.1, 0.2, 0.7], [0.8, 0.05, 0.15]], dtype=torch.float)
        real_idxes = torch.tensor([2, 2], dtype=torch.long)
        tac = TopkAccCalculator()
        tac.extend_pred_and_real(pred_probs_list, real_idxes)
        self.assertEqual(tac.get_topk_acc(1), 0.5)

    def test_top2_acc(self):
        pred_probs_list = torch.tensor([[0.1, 0.2, 0.7], [0.8, 0.05, 0.15]], dtype=torch.float)
        real_idxes = torch.tensor([2, 2], dtype=torch.long)
        tac = TopkAccCalculator()
        tac.extend_pred_and_real(pred_probs_list, real_idxes)
        self.assertEqual(tac.get_topk_acc(2), 1.0)


if __name__ == "__main__":
    pass
