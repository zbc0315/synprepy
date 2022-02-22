#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Time    : 2022/1/11 17:23
# @Author  : zhangbc0315@outlook.com
# @File    : model_utils.py
# @Software: PyCharm

import os
import logging

import torch


class ModelUtils:

    @classmethod
    def _get_model_fn(cls, epoch: int, train_loss: float, train_acc: float, test_loss: float, test_acc: float):
        return f"model_e{epoch}_al{train_loss}_aa{train_acc}_el{test_loss}_ea{test_acc}_.pt"

    @classmethod
    def save_model(cls, dp: str, model: torch.nn.Module, epoch: int,
                   train_loss: float, train_acc: float, test_loss: float, test_acc: float):
        model_fn = cls._get_model_fn(epoch, train_loss, train_acc, test_loss, test_acc)
        torch.save({"model": model.state_dict()}, os.path.join(dp, model_fn))

    @classmethod
    def _get_test_loss(cls, fn: str) -> float:
        for part in fn.split('_'):
            if part.startswith('ea'):
                return float(part[2:])
        raise ValueError("模型命名格式有错")

    @classmethod
    def load_best_model(cls, dp: str):
        sorted_fns = sorted(os.listdir(dp), key=cls._get_test_loss)
        fp = os.path.join(dp, sorted_fns[-1])
        logging.info(f"load model: {fp}")
        return torch.load(fp)["model"]


if __name__ == "__main__":
    pass
