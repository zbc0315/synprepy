#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Time    : 2021/8/25 14:13
# @Author  : zhangbc0315@outlook.com
# @File    : gcn_model.py
# @Software: PyCharm

import os

import torch
import torch.nn as nn
import torch.nn.functional as F
import torch_geometric.nn as gnn

from nn_utils.gcn_base_model import NodeModel, EdgeModel, GlobalModel
from config import RDPath


class GcnModel(nn.Module):

    def __init__(self, num_node_features, num_edge_features, num_global_features):
        super(GcnModel, self).__init__()
        self.num_global_features = num_global_features
        self.meta_r1 = gnn.MetaLayer(EdgeModel(num_node_features, num_edge_features, 256),
                                     NodeModel(num_node_features, 256, 128),
                                     GlobalModel(128, num_global_features, 128))
        self.meta_r2 = gnn.MetaLayer(EdgeModel(128, 256, 256),
                                     NodeModel(128, 256, 128),
                                     GlobalModel(128, 128, 128))
        self.meta_r3 = gnn.MetaLayer(EdgeModel(128, 256, 256),
                                     NodeModel(128, 256, 128),
                                     GlobalModel(128, 128, 128))
        self.meta_r4 = gnn.MetaLayer(EdgeModel(128, 256, 256),
                                     NodeModel(128, 256, 128),
                                     GlobalModel(128, 128, 128))

        self.meta_p1 = gnn.MetaLayer(EdgeModel(num_node_features, num_edge_features, 256),
                                     NodeModel(num_node_features, 256, 128),
                                     GlobalModel(128, 0, 128))
        self.meta_p2 = gnn.MetaLayer(EdgeModel(128, 256, 256),
                                     NodeModel(128, 256, 128),
                                     GlobalModel(128, 128, 128))
        self.meta_p3 = gnn.MetaLayer(EdgeModel(128, 256, 256),
                                     NodeModel(128, 256, 128),
                                     GlobalModel(128, 128, 128))
        self.meta_p4 = gnn.MetaLayer(EdgeModel(128, 256, 256),
                                     NodeModel(128, 256, 128),
                                     GlobalModel(128, 128, 128))

        self.r_lin1 = nn.Linear(128, 128)
        # self.r_lin2 = nn.Linear(128, 128)
        # self.r_lin3 = nn.Linear(128, 128)

        self.p_lin1 = nn.Linear(128, 128)
        # self.p_lin2 = nn.Linear(128, 128)
        # self.p_lin3 = nn.Linear(128, 128)

    def forward(self, r_data, p_data):
        rx, re, rc, rb, rg = r_data.x, r_data.edge_attr, r_data.edge_index, r_data.batch, r_data.g
        px, pe, pc, pb, pg = p_data.x, p_data.edge_attr, p_data.edge_index, p_data.batch, None

        if self.num_global_features == 0:
            rg = None
        rx, re, rg = self.meta_r1(rx, rc, re, rg, rb)
        rx, re, rg = self.meta_r2(rx, rc, re, rg, rb)
        rx, re, rg = self.meta_r3(rx, rc, re, rg, rb)
        rx, re, rg = self.meta_r4(rx, rc, re, rg, rb)

        px, pe, pg = self.meta_p1(px, pc, pe, pg, pb)
        px, pe, pg = self.meta_p2(px, pc, pe, pg, pb)
        px, pe, pg = self.meta_p3(px, pc, pe, pg, pb)
        px, pe, pg = self.meta_p4(px, pc, pe, pg, pb)

        # rg = F.relu(self.r_lin1(rg))
        # rg = F.relu(self.r_lin2(rg))
        rg = self.r_lin1(rg)

        # pg = F.relu(self.p_lin1(pg))
        # pg = F.relu(self.p_lin2(pg))
        pg = self.p_lin1(pg)

        return F.sigmoid(F.cosine_similarity(rg, pg))

    @classmethod
    def get_model_and_opt_state(cls, with_spec: bool):
        if with_spec:
            dp = RDPath.MODEL_CAT_SOL_DP
        else:
            dp = RDPath.MODEL_NO_CAT_SOL_DP

        fns = list(os.listdir(dp))
        fns = sorted(fns, key=lambda x: float(x.split('_')[3][1:]))
        best_fn = fns[-1]
        print(f'load model from {dp} : {best_fn}')
        state = torch.load(os.path.join(dp, best_fn))
        return state['model'], state['opt']


if __name__ == "__main__":
    pass
