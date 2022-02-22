# -*- coding: utf-8 -*-
# @Time     : 2021/5/17 21:33
# @Author   : zbc
# @Email    : zhangbc0315@outlook.com
# @FileName : eval_and_predict.py
# @Software : PyCharm

import os

import pandas as pd
import torch

# from log_utils import Logger
from config import SSPath
from rxn_spectators_selector.ss_model import SSModel
from rxn_spectators_selector.rxn_data import RxnData


class EvalAndPredict:

    _device = "cuda"
    _spec_count_filter_fps = {"cat": SSPath.CATS_COUNT_FILTER_FP,
                              "sol": SSPath.SOLS_COUNT_FILTER_FP}
    _spec_topk_fps = {"cat": SSPath.CAT_GNN_TOPK_FP,
                      "sol": SSPath.SOL_GNN_TOPK_FP}
    _spec_topk_acc_fps = {"cat": SSPath.CAT_GNN_TOPK_ACC_FP,
                          "sol": SSPath.SOL_GNN_TOPK_ACC_FP}
    _spec_topk_fp = None
    _spec_topk_acc_fp = None
    _spec_count_filter_fp = None
    _spec_count_filter_df = None
    _model = None

    @classmethod
    def _init(cls, model_type: str):
        # Logger.info("Eval And Predict Initializing")
        cls._spec_topk_fp = cls._spec_topk_fps[model_type]
        cls._spec_topk_acc_fp = cls._spec_topk_acc_fps[model_type]
        cls._spec_count_filter_fp = cls._spec_count_filter_fps[model_type]
        cls._spec_count_filter_df = pd.read_csv(cls._spec_count_filter_fp, sep='\t', encoding="utf-8")
        cls._model = SSModel.load_model(model_type).to(cls._device)
        cls._model.eval()
        # Logger.info("Eval And Predict Initialized")

    @classmethod
    def _init_if_necessary(cls, model_type: str):
        if cls._spec_count_filter_df is None:
            cls._init(model_type)

    @classmethod
    def get_sid_from_oid(cls, oid: int) -> int:
        """

        :param oid: one-hot id
        :return: spectator id
        """
        return cls._spec_count_filter_df.loc[oid, "sid"]

    @classmethod
    def _calc_eval(cls, predict, real):
        sort_pred_idxes = torch.sort(predict, dim=1, descending=True)[1]
        max_real_idxes = real.max(dim=1)[1]
        res = []
        for sort_pred_idx, max_real_idx in zip(sort_pred_idxes, max_real_idxes):
            top_k = list(sort_pred_idx).index(max_real_idx)
            res.append(top_k)
        return res

    @classmethod
    def eval_one_batch(cls, r_batch_data, p_batch_data):
        r_batch_data = r_batch_data.to(cls._device)
        p_batch_data = p_batch_data.to(cls._device)
        predict = cls._model(r_batch_data, p_batch_data)
        return cls._calc_eval(predict, r_batch_data.y)

    @classmethod
    def _write_topk(cls, model_type: str):
        data_set = RxnData(model_type, "test")
        top_ks = []
        for n, (r_batch_data, p_batch_data) in enumerate(data_set.get_batch_data()):
            # if n % 10 == 0:
            #     Logger.info(f"eval: {n}")
            top_ks.extend(cls.eval_one_batch(r_batch_data, p_batch_data))
        top_ks = [str(i) for i in top_ks]
        with open(cls._spec_topk_fp, 'w', encoding="utf-8")as f:
            f.write('\n'.join(top_ks))

    @classmethod
    def _write_topk_acc(cls, res):
        with open(cls._spec_topk_acc_fp, 'w', encoding='utf-8')as f:
            f.write("k\ttopk accuracy\n")
            for k, acc in res.items():
                f.write(f"{k}\t{acc}\n")

    @classmethod
    def _analyze_topk(cls, max_k):
        with open(cls._spec_topk_fp, 'r', encoding="utf-8")as f:
            top_ks = [int(k) for k in f.read().split('\n')]
        top_k_count = {}
        current_count = 0
        for i in range(max_k):
            if i in top_ks:
                current_count += top_ks.count(i)
            top_k_count[i+1] = current_count / len(top_ks)
            print(f"top {i+1}: {top_k_count[i+1]}")
        cls._write_topk_acc(top_k_count)

    @classmethod
    def eval(cls, model_type: str):
        cls._init(model_type)
        # if not os.path.exists(cls._spec_topk_fp):
        cls._write_topk(model_type)
        cls._analyze_topk(10)


if __name__ == '__main__':
    EvalAndPredict.eval("sol")
