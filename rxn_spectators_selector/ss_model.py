# -*- coding: utf-8 -*-
# @Time     : 2021/4/28 16:01
# @Author   : zbc
# @Email    : zhangbc0315@outlook.com
# @FileName : ss_model.py
# @Software : PyCharm

import os

import torch
import torch.nn as nn
from torch.nn import Module
import torch.nn.functional as F
import torch_geometric.nn as gnn

from nn_utils.gcn_base_model import NodeModel, EdgeModel, GlobalModel
from config.path import SpectatorsSelectorPath as SSPath


class SSModel(Module):

    _model_dps = {"cat": SSPath.CAT_GNN_MODEL_DP,
                  "sol": SSPath.SOL_GNN_MODEL_DP}

    @classmethod
    def _get_correct_value_from_fn(cls, fn: str) -> float:
        return float(fn[:-3].split('_')[-1][1:])

    @classmethod
    def _get_best_mode_fp(cls, model_dp: str) -> str:
        fns = list(os.listdir(model_dp))
        fn_to_c = {fn: cls._get_correct_value_from_fn(fn) for fn in fns}
        sorted_fns = sorted(fn_to_c.keys(), key=lambda x: fn_to_c[x], reverse=False)
        print(sorted_fns[-1])
        return os.path.join(model_dp, sorted_fns[-1])

    @classmethod
    def load_model(cls, model_type: str) -> Module:
        model_fp = cls._get_best_mode_fp(cls._model_dps[model_type])
        return torch.load(model_fp)
    
    def __init__(self, num_node_features, num_edge_features, num_global, out_channels):
        super(SSModel, self).__init__()
        self.out_channels = out_channels
        # num_global = 0
        self.meta_r1 = gnn.MetaLayer(EdgeModel(num_node_features, num_edge_features, 512),
                                     NodeModel(num_node_features, 512, 128),
                                     GlobalModel(128, num_global, 128))
        self.meta_r2 = gnn.MetaLayer(EdgeModel(128, 512, 512),
                                     NodeModel(128, 512, 128),
                                     GlobalModel(128, 128, 128))
        self.meta_r3 = gnn.MetaLayer(EdgeModel(128, 512, 512),
                                     NodeModel(128, 512, 128),
                                     GlobalModel(128, 128, 128))
        self.meta_r4 = gnn.MetaLayer(EdgeModel(128, 512, 128),
                                     NodeModel(128, 128, 128),
                                     GlobalModel(128, 128, 256))

        self.meta_p1 = gnn.MetaLayer(EdgeModel(num_node_features, num_edge_features, 512),
                                     NodeModel(num_node_features, 512, 128),
                                     GlobalModel(128, num_global, 128))
        self.meta_p2 = gnn.MetaLayer(EdgeModel(128, 512, 512),
                                     NodeModel(128, 512, 128),
                                     GlobalModel(128, 128, 128))
        self.meta_p3 = gnn.MetaLayer(EdgeModel(128, 512, 512),
                                     NodeModel(128, 512, 128),
                                     GlobalModel(128, 128, 128))
        self.meta_p4 = gnn.MetaLayer(EdgeModel(128, 512, 128),
                                     NodeModel(128, 128, 128),
                                     GlobalModel(128, 128, 256))

        self.lin1 = nn.Linear(512, 512)
        self.lin2 = nn.Linear(512, out_channels)

    def forward(self, r_batch_data, p_batch_data):
        rx, re, rc, rb, rg = r_batch_data.x, r_batch_data.edge_attr, r_batch_data.edge_index, r_batch_data.batch, r_batch_data.g
        px, pe, pc, pb, pg = p_batch_data.x, p_batch_data.edge_attr, p_batch_data.edge_index, p_batch_data.batch, p_batch_data.g
        # rx, re, rc, px, pe, pc, batch = data.rx, data.re, data.rc, data.px, data.pe, data.pc, data.batch
        rx, re, rg = self.meta_r1(rx, rc, re, rg, rb)
        rx, re, rg = self.meta_r2(rx, rc, re, rg, rb)
        rx, re, rg = self.meta_r3(rx, rc, re, rg, rb)
        rx, re, rg = self.meta_r4(rx, rc, re, rg, rb)

        px, pe, pg = self.meta_p1(px, pc, pe, pg, pb)
        px, pe, pg = self.meta_p2(px, pc, pe, pg, pb)
        px, pe, pg = self.meta_p3(px, pc, pe, pg, pb)
        px, pe, pg = self.meta_p4(px, pc, pe, pg, pb)

        g = torch.cat([rg, pg], dim=1)
        y = F.relu(self.lin1(g))
        return F.log_softmax(self.lin2(y), dim=-1)


if __name__ == '__main__':
    pass
