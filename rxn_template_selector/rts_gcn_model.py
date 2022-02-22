#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Time    : 2021/8/10 12:37
# @Author  : zhangbc0315@outlook.com
# @File    : rts_gcn_model.py
# @Software: PyCharm

import os

import torch
import torch.nn as nn
import torch.nn.functional as F
import torch_geometric.nn as gnn
import torch_geometric.data as gdata
from torch_geometric.nn import BatchNorm

from nn_utils.gcn_base_model import NodeModel, EdgeModel, GlobalModel
from nn_utils import ModelNameUtils
from config import BTPath


class RtsGcnModel(nn.Module):

    cbe_tids = None
    cbe_loss_tids = None

    def __init__(self, num_node_features: int, num_edge_features: int, num_out_channels: int):
        super(RtsGcnModel, self).__init__()
        print(f"node: {num_node_features}, edge: {num_edge_features}, out: {num_out_channels}")
        num_node_features -= 2
        # num_edge_features -= 1
        self.num_out_channels = num_out_channels
        # num_out_channels = self._get_num_y()
        self.node_normal = gnn.BatchNorm(num_node_features)
        self.edge_normal = gnn.BatchNorm(num_edge_features)
        self.meta1 = gnn.MetaLayer(EdgeModel(num_node_features, num_edge_features, 512),
                                   NodeModel(num_node_features, 512, 128),
                                   GlobalModel(128, 0, 128))
        self.meta2 = gnn.MetaLayer(EdgeModel(128, 512, 512),
                                   NodeModel(128, 512, 128),
                                   GlobalModel(128, 128, 128))
        self.meta3 = gnn.MetaLayer(EdgeModel(128, 512, 512),
                                   NodeModel(128, 512, 128),
                                   GlobalModel(128, 128, 128))
        self.meta4 = gnn.MetaLayer(EdgeModel(128, 512, 512),
                                   NodeModel(128, 512, 128),
                                   GlobalModel(128, 128, 128))
        self.meta5 = gnn.MetaLayer(EdgeModel(128, 512, 512),
                                   NodeModel(128, 512, 128),
                                   GlobalModel(128, 128, 128))
        self.meta6 = gnn.MetaLayer(EdgeModel(128, 512, 128),
                                   NodeModel(128, 128, 128),
                                   GlobalModel(128, 128, 256))
        self.lin1 = nn.Linear(256, 512)
        self.lin2 = nn.Linear(512, num_out_channels)

    def forward(self, data: gdata):
        x, edge_index, e, batch = data.x, data.edge_index, data.edge_attr, data.batch
        # e[:, 0] = e[:, -1]
        # e = e[:, :-1]
        x = x[:, :-2]
        x = self.node_normal(x)
        e = self.edge_normal(e)
        x, e, g = self.meta1(x, edge_index, e, None, batch)
        x, e, g = self.meta2(x, edge_index, e, g, batch)
        x, e, g = self.meta3(x, edge_index, e, g, batch)
        # x, e, g = self.meta4(x, edge_index, e, g, batch)
        x, e, g = self.meta5(x, edge_index, e, g, batch)
        x, e, g = self.meta6(x, edge_index, e, g, batch)

        y = F.elu(self.lin1(g))
        return g, self.lin2(y)

    @classmethod
    def get_model_and_opt(cls, model_dp: str):
        model_fns = os.listdir(model_dp)
        if len(model_fns) == 0:
            print(f"没有模型可供加载")
            return None, None, 0
        model_fns = sorted(model_fns, key=lambda x: ModelNameUtils.parse_model_name(x)['loss'])
        model_fn = model_fns[0]
        model_fp = os.path.join(model_dp, model_fn)
        res = torch.load(model_fp)
        print(f"加载模型：{model_fn}")
        return res['model'], res['opt'], ModelNameUtils.parse_model_name(model_fn)['idx']

    @classmethod
    def get_model_fp(cls, models_dp: str, sort_key: str = "best"):
        return "C:\\Users\\zhang\\Documents\\Data\\RxnPredictor\\GNN\\BaseTemplate_AtomicNum\\model_neg\\model_54_0.020973628762365443_0.4128210837928821.pt"
        """ 在highway gcn中，作为预训练gcn描述符提供者时使用

        :param models_dp:
        :param sort_key: "best", "latest"
        :return:
        """
        key_to_idx = {"best": 3,
                      "latest": 1}
        model_fns = os.listdir(models_dp)
        if len(model_fns) == 0:
            raise AttributeError(f"加载模型：{models_dp} 中不存在任何文件")
        model_fns = sorted(model_fns, key=lambda x: float(x.split('.pt')[0].split('_')[key_to_idx[sort_key]]))
        latest_model_fn = model_fns[0]
        print(f"model: {latest_model_fn}")
        return os.path.join(models_dp, latest_model_fn)

    def _load_cbe_tids(self):
        if self.cbe_tids is None:
            with open(BTPath.CBE_TIDS_FP, 'r', encoding='utf-8')as f:
                self.cbe_tids = [int(tid) for tid in f.read().split('\n')]
                self.cbe_loss_tids = set(range(self.num_out_channels)) - set(self.cbe_tids)

    def _get_num_y(self):
        self._load_cbe_tids()
        return len(self.cbe_tids)

    def change_y(self, y):
        return y
        self._load_cbe_tids()
        return y[:, self.cbe_tids]


if __name__ == "__main__":
    pass
