#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Time    : 2021/8/2 14:33
# @Author  : zhangbc0315@outlook.com
# @File    : model_gcn.py
# @Software: PyCharm

import torch
from torch import nn
import torch.nn.functional as F
from torch.nn import Sequential as Seq, Linear as Lin, ReLU as ReLU
from torch_scatter import scatter_mean, scatter_sum
import torch_geometric.nn as gnn

from config import BTPath


class OldBatchNorm(nn.BatchNorm1d):

    """ copy from torch_geometric/nn/norm/batch_norm.py 1.6.0

    """

    def __init__(self, in_channels, eps=1e-5, momentum=0.1, affine=True,
                 track_running_stats=True):
        super(OldBatchNorm, self).__init__(in_channels, eps, momentum, affine,
                                           track_running_stats)

    def forward(self, x):
        """"""
        return super(OldBatchNorm, self).forward(x)

    def __repr__(self):
        return ('{}({}, eps={}, momentum={}, affine={}, '
                'track_running_stats={})').format(self.__class__.__name__,
                                                  self.num_features, self.eps,
                                                  self.momentum, self.affine,
                                                  self.track_running_stats)


class EdgeModel(nn.Module):
    def __init__(self, num_node_features, num_edge_features, out_features):
        super(EdgeModel, self).__init__()
        self.edge_mlp = Seq(Lin(num_node_features+num_node_features+num_edge_features, 128), ReLU(), Lin(128, 128), ReLU(), Lin(128, out_features))

    def forward(self, src, dest, edge_attr, u, batch):
        # source, target: [E, F_x], where E is the number of edges.
        # edge_attr: [E, F_e]
        # u: [B, F_u], where B is the number of graphs.
        # batch: [E] with max entry B - 1.
        out = torch.cat([src, dest, edge_attr], 1)
        return self.edge_mlp(out)


class NodeModel(torch.nn.Module):
    def __init__(self, num_node_features, num_edge_features_out, out_features):
        super(NodeModel, self).__init__()
        self.node_mlp_1 = Seq(Lin(num_node_features + num_edge_features_out, 256), ReLU(), Lin(256, 256), ReLU(), Lin(256, 256))
        self.node_mlp_2 = Seq(Lin(num_node_features + 256, 256), ReLU(), Lin(256, out_features))

    def forward(self, x, edge_index, edge_attr, u, batch):
        # x: [N, F_x], where N is the number of nodes.
        # edge_index: [2, E] with max entry N - 1.
        # edge_attr: [E, F_e]
        # u: [B, F_u]
        # batch: [N] with max entry B - 1.
        row, col = edge_index
        out = torch.cat([x[row], edge_attr], dim=1)
        out = self.node_mlp_1(out)
        out = scatter_mean(out, col, dim=0, dim_size=x.size(0))
        out = torch.cat([x, out], dim=1)
        return self.node_mlp_2(out)


class GlobalModel(torch.nn.Module):
    def __init__(self, num_node_features, num_global_features, out_channels):
        super(GlobalModel, self).__init__()
        self.global_mlp = Seq(Lin(num_global_features+num_node_features, 256), ReLU(), Lin(256, 256), ReLU(), Lin(256, out_channels))

    def forward(self, x, edge_index, edge_attr, u, batch):
        # x: [N, F_x], where N is the number of nodes.
        # edge_index: [2, E] with max entry N - 1.
        # edge_attr: [E, F_e]
        # u: [B, F_u]
        # batch: [N] with max entry B - 1.
        out = torch.cat([u, scatter_mean(x, batch, dim=0)], dim=1)
        return self.global_mlp(out)


class Net(nn.Module):
    def __init__(self, num_node_features: int, num_edge_features: int, device: str):
        super(Net, self).__init__()
        self.bn_node = OldBatchNorm(num_node_features)
        self.bn_edge = OldBatchNorm(num_edge_features)
        self.device = device
        self.meta1 = gnn.MetaLayer(EdgeModel(num_node_features, num_edge_features, 512),
                                   NodeModel(num_node_features, 512, 128),
                                   GlobalModel(128, 1, 128))
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
                                   None,
                                   None)
        self.lin1 = nn.Linear(128, 128)
        self.lin2 = nn.Linear(128, 1)

    def forward(self, data, batch_size):
        x, edge_index, e, g, edge_type = data.x, data.edge_index, data.edge_attr, data.g, data.edge_type
        x = self.bn_node(x)
        e = self.bn_edge(e)
        x, e, g = self.meta1(x, edge_index, e, g, data.batch)
        x, e, g = self.meta2(x, edge_index, e, g, data.batch)
        x, e, g = self.meta3(x, edge_index, e, g, data.batch)
        x, e, g = self.meta4(x, edge_index, e, g, data.batch)
        x, e, g = self.meta5(x, edge_index, e, g, data.batch)
        _, e, _ = self.meta6(x, edge_index, e, g, data.batch)

        eb = self.get_edge_batch_tensor(batch_size, self.device)

        selected = e.t()[:, edge_type == 0].t()
        selected = scatter_sum(selected, eb, dim=0)
        y = selected
        y = self.lin1(y)
        y = self.lin2(F.relu(y))
        return y

    @classmethod
    def get_edge_batch_tensor(cls, batch_size, device):
        arr = []
        for i in range(batch_size):
            arr.append(i)
            arr.append(i)
        return torch.tensor(arr, dtype=torch.long, device=device)

    @classmethod
    def get_model(cls, device: str):
        model = cls(11, 2, device).to(device)
        state_dict = torch.load(BTPath.BOND_ENERGY_MODEL)['model']
        model.load_state_dict(torch.load(BTPath.BOND_ENERGY_MODEL)['model'])
        return model


if __name__ == "__main__":
    pass
