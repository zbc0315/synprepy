# -*- coding: utf-8 -*-
# @Time     : 2021/4/28 15:59
# @Author   : zbc
# @Email    : zhangbc0315@outlook.com
# @FileName : gcn_base_model.py
# @Software : PyCharm

import torch
import torch_scatter
from torch.nn import Module, Sequential, Linear, ReLU, ELU, Dropout


class EdgeModel(Module):
    def __init__(self, num_node_features, num_edge_features, out_features):
        super(EdgeModel, self).__init__()
        self.edge_mlp = Sequential(Linear(num_node_features + num_node_features + num_edge_features, 128),
                                   ELU(),
                                   # Dropout(0.3),
                                   Linear(128, out_features))

    def forward(self, src, dest, edge_attr, u, batch):
        out = torch.cat([src, dest, edge_attr], 1)
        return self.edge_mlp(out)


class NodeModel(Module):
    def __init__(self, num_node_features, num_edge_features_out, out_features):
        super(NodeModel, self).__init__()
        self.node_mlp_1 = Sequential(Linear(num_node_features + num_edge_features_out, 256),
                                     ELU(),
                                     # Dropout(0.3),
                                     Linear(256, 256))
        self.node_mlp_2 = Sequential(Linear(num_node_features + 256, 256),
                                     ELU(),
                                     # Dropout(0.3),
                                     Linear(256, out_features))

    def forward(self, x, edge_index, edge_attr, u, batch):
        row, col = edge_index
        out = torch.cat([x[row], edge_attr], dim=1)
        out = self.node_mlp_1(out)
        out = torch_scatter.scatter_mean(out, col, dim=0, dim_size=x.size(0))
        out = torch.cat([x, out], dim=1)
        return self.node_mlp_2(out)


class GlobalModel(Module):
    def __init__(self, num_node_features, num_global_features, out_channels):
        super(GlobalModel, self).__init__()
        self.global_mlp = Sequential(Linear(num_global_features + num_node_features, 256),
                                     ELU(),
                                     # Dropout(0.3),
                                     Linear(256, out_channels))

    def forward(self, x, edge_index, edge_attr, u, batch):
        if u is None:
            out = torch_scatter.scatter_mean(x, batch, dim=0)
        else:
            out = torch.cat([u, torch_scatter.scatter_mean(x, batch, dim=0)], dim=1)
            # print(f"size: {out.size()}")
        return self.global_mlp(out)


if __name__ == '__main__':
    pass
