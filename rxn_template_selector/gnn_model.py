#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Time    : 2022/1/9 17:47
# @Author  : zhangbc0315@outlook.com
# @File    : yp_gnn_model.py
# @Software: PyCharm

import logging

import torch
from torch import nn
import torch.nn.functional as F
import torch_geometric.nn as gnn

from utils.gnn_base_model import NodeModel, EdgeModel, GlobalModel


class GNNModel(nn.Module):
    
    def __init__(self, num_node_features, num_edge_features, num_global_features, out_channels):
        super(GNNModel, self).__init__()
        logging.info(f"num node features: {num_node_features}, "
                     f"num edge features: {num_edge_features}, "
                     f"num global features: {num_global_features}, "
                     f"num out channels: {out_channels}")
        # num_node_features -= 2
        # act_func = nn.ELU
        self.node_normal = gnn.BatchNorm(num_node_features)
        self.edge_normal = gnn.BatchNorm(num_edge_features)
        # self.node_layer = nn.Linear(num_node_features, num_node_features)
        # self.edge_layer = nn.Linear(num_edge_features, num_edge_features)
        self.meta1 = gnn.MetaLayer(EdgeModel(num_node_features, num_edge_features, 512),
                                   NodeModel(num_node_features, 512, 128),
                                   GlobalModel(128, num_global_features, 128))
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
        self.lin2 = nn.Linear(512, out_channels)

    def forward(self, data):
        # x, edge_index, e, g, batch = data.x, data.edge_index, data.edge_attr, data.g, data.batch
        x, edge_index, e, g, batch = data.x, data.edge_index, data.edge_attr, data.g, data.batch
        # x = x[:, :-2]
        x = self.node_normal(x)
        e = self.edge_normal(e)

        x, e, g = self.meta1(x, edge_index, e, g, batch)
        x = F.dropout(x, 0.1)
        x, e, g = self.meta2(x, edge_index, e, g, batch)
        x = F.dropout(x, 0.1)
        x, e, g = self.meta3(x, edge_index, e, g, batch)
        x = F.dropout(x, 0.1)
        x, e, g = self.meta4(x, edge_index, e, g, batch)
        x = F.dropout(x, 0.1)
        x, e, g = self.meta5(x, edge_index, e, g, batch)
        x = F.dropout(x, 0.1)
        x, e, g = self.meta6(x, edge_index, e, g, batch)

        y = F.elu(self.lin1(g))
        return g, self.lin2(y)


if __name__ == "__main__":
    pass
