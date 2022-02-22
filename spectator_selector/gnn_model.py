#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Time    : 2022/2/22 10:45
# @Author  : zhangbc0315@outlook.com
# @File    : gnn_model.py
# @Software: PyCharm

import torch
from torch import nn
import torch.nn.functional as F
import torch_geometric.nn as gnn

from utils.gnn_base_model import NodeModel, EdgeModel, GlobalModel


class CSSSGNNModel(nn.Module):

    def __init__(self, num_node_features, num_edge_features, num_global_features):
        super(CSSSGNNModel, self).__init__()
        self.node_normal = gnn.BatchNorm(num_node_features)
        self.edge_normal = gnn.BatchNorm(num_edge_features)
        # act_func = nn.LeakyReLU
        self.meta_r1 = gnn.MetaLayer(EdgeModel(num_node_features, num_edge_features, 512),
                                     NodeModel(num_node_features, 512, 128),
                                     GlobalModel(128, num_global_features, 128))
        self.meta_r2 = gnn.MetaLayer(EdgeModel(128, 512, 512),
                                     NodeModel(128, 512, 128),
                                     GlobalModel(128, 128, 128))
        self.meta_r3 = gnn.MetaLayer(EdgeModel(128, 512, 512),
                                     NodeModel(128, 512, 128),
                                     GlobalModel(128, 128, 128))
        self.meta_r4 = gnn.MetaLayer(EdgeModel(128, 512, 512),
                                     NodeModel(128, 512, 128),
                                     GlobalModel(128, 128, 128))
        self.meta_r5 = gnn.MetaLayer(EdgeModel(128, 512, 512),
                                     NodeModel(128, 512, 128),
                                     GlobalModel(128, 128, 256))

        self.meta_p1 = gnn.MetaLayer(EdgeModel(num_node_features, num_edge_features, 512),
                                     NodeModel(num_node_features, 512, 128),
                                     GlobalModel(128, 0, 128))
        self.meta_p2 = gnn.MetaLayer(EdgeModel(128, 512, 512),
                                     NodeModel(128, 512, 128),
                                     GlobalModel(128, 128, 128))
        self.meta_p3 = gnn.MetaLayer(EdgeModel(128, 512, 512),
                                     NodeModel(128, 512, 128),
                                     GlobalModel(128, 128, 256))
        # self.meta_p4 = gnn.MetaLayer(EdgeModel(128, 512, 512),
        #                              NodeModel(128, 512, 128),
        #                              GlobalModel(128, 128, 128))
        self.r_lin1 = nn.Linear(256, 256)
        self.p_lin1 = nn.Linear(256, 256)
        self.y_lin = nn.Linear(512, 1)

    def forward(self, r_data, p_data):
        rx, re, rc, rb, rg = r_data.x, r_data.edge_attr, r_data.edge_index, r_data.batch, r_data.g
        px, pe, pc, pb, pg = p_data.x, p_data.edge_attr, p_data.edge_index, p_data.batch, None

        rx = self.node_normal(rx)
        px = self.node_normal(px)
        re = self.edge_normal(re)
        pe = self.edge_normal(pe)
        rx, re, rg = self.meta_r1(rx, rc, re, rg, rb)
        # rx = F.dropout(rx, 0.1)
        rx, re, rg = self.meta_r2(rx, rc, re, rg, rb)
        # rx = F.dropout(rx, 0.1)
        rx, re, rg = self.meta_r3(rx, rc, re, rg, rb)
        # rx = F.dropout(rx, 0.1)
        rx, re, rg = self.meta_r4(rx, rc, re, rg, rb)
        # rx = F.dropout(rx, 0.1)
        rx, re, rg = self.meta_r5(rx, rc, re, rg, rb)

        px, pe, pg = self.meta_p1(px, pc, pe, pg, pb)
        px, pe, pg = self.meta_p2(px, pc, pe, pg, pb)
        px, pe, pg = self.meta_p3(px, pc, pe, pg, pb)

        rg = F.elu(self.r_lin1(rg))
        pg = F.elu(self.p_lin1(pg))
        g = torch.cat([rg, pg], 1)
        return F.sigmoid(self.y_lin(g)) * 100


if __name__ == "__main__":
    pass
