#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Time    : 2021/11/10 15:15
# @Author  : zhangbc0315@outlook.com
# @File    : predict.py
# @Software: PyCharm


import torch
from torch_geometric.data import Batch

from rxn_determiner.gcn_model import GcnModel
from nn_utils.gnn_data_utils import GnnDataUtils


class Predict:
    cats_len = 325
    sols_len = 321
    device = 'cpu'

    no_spec = 'model_opt_89_a0.9439160950969608_p0.9386040860236788_r0.8961484399375139_.pt'
    # no_spec = 'model_opt_32_a0.9490094499449064_p0.938478012435774_r0.8946966016944882_.pt'
    # no_spec = 'model_opt_51_a0.9556703671391639_p0.946675693701967_r0.9086443344238879_.pt'

    with_spec = 'model_opt_74_a0.9702884915741274_p0.9720672494896699_r0.9409673341646233_.pt'

    _model_fps = {'no_spectators': f'C:\\Users\\zhang\\Documents\\Data\\RxnMolPredictor\\RxnDeterminer\\model_no_cat_sol\\{no_spec}',
                  'with_spectators': f'C:\\Users\\zhang\\Documents\\Data\\RxnMolPredictor\\RxnDeterminer\\model_cat_sol\\{with_spec}'}

    @classmethod
    def load_model(cls, with_spectators: bool):
        if with_spectators:
            num_global_feature = 650
            model_fp = cls._model_fps['with_spectators']
        else:
            num_global_feature = 0
            model_fp = cls._model_fps['no_spectators']
        model = GcnModel(6, 12, num_global_feature)
        model.load_state_dict(torch.load(model_fp)['model'])
        return model

    @classmethod
    def get_global_feature(cls, cat_idx: int, sol_idx: int):
        cat_vector = [0]*(cls.cats_len + 2)
        cat_vector[cat_idx + 1] = 1
        sol_vector = [0]*(cls.sols_len + 2)
        sol_vector[sol_idx + 1] = 1
        cat_vector.extend(sol_vector)
        return [cat_vector]

    @classmethod
    def get_data(cls, rs_smiles: [str], p_smiles: str, cat_idx: int, sol_idx: int):
        rs_data = GnnDataUtils.gnn_data_from_smileses(rs_smiles)
        p_data = GnnDataUtils.gnn_data_from_smiles(p_smiles)
        global_feature = cls.get_global_feature(cat_idx, sol_idx)
        rs_data.g = torch.tensor(torch.tensor(global_feature), dtype=torch.float)
        return Batch.from_data_list([rs_data]), Batch.from_data_list([p_data])

    @classmethod
    def predict_by_smiles(cls, with_spectators: bool, rs_smiles: [str], p_smiles: str, cat_idx: int, sol_idx: int):
        model = cls.load_model(with_spectators)
        rs_data, p_data = cls.get_data(rs_smiles, p_smiles, cat_idx, sol_idx)
        rs_data = rs_data.to(cls.device)
        p_data = p_data.to(cls.device)
        pred = model(rs_data, p_data)
        print(pred)


if __name__ == "__main__":
    rs = ['CC1=CC=CN=C1C2=CC=C(OC)C=C2', 'ClN1C(CCC1=O)=O']
    p1 = 'ClC1=C(OC)C=CC(C2=NC=CC=C2C)=C1'
    p2 = 'ClC(C=C(OC)C=C1)=C1C2=NC=CC=C2C'
    rs = ['C1(N2CCCCC2)=CCCCC1', 'O=NC1=CC=CC=C1']
    p1 = 'O=C1C(N(O)C2=CC=CC=C2)CCCC1'
    p2 = 'O=C1C(ONC2=CC=CC=C2)CCCC1'

    rs = ['O=C(C=C1)N(C2=CC=CC=C2)C1=O', 'BrCC(C(OCC)=O)=CC']
    p1 = 'O=C1N(C2=CC=CC=C2)C(C3C=C(C(OCC)=O)C(C)C31)=O'
    p2 = 'O=C1N(C2=CC=CC=C2)C(C3CC=C(C(OCC)=O)CC31)=O'
    c1 = 98
    c2 = 89
    l1 = l2 = 1

    rs = ['O=C(N)OCC/C=C/CC']
    p1 = 'O=C1OCCC2N1C2CC'
    p2 = 'O=C1OCC(/C=C/CC)N1'
    c1 = 55
    c2 = 325
    l1 = l2 = 300
    #
    # rs = ['O=C/C=C/C1=CC=CC=C1', 'O=C(C(C)(C)O)/C=C/C(OCC)=O']
    # p1 = 'OC(C)(C)C1(O2)CC(C(OCC)=O)C(C3=CC=CC=C3)C1C2=O'
    # p2 = 'OC(C)(C)C(O1)=CC(C(OCC)=O)C(CC2=CC=CC=C2)C1=O'
    # c1 = 188
    # c2 = 64
    # l1 = l2 = 314

    # rs = ['O=S(C1=C2C=CC=C1)(N=C2C3=CC=CC=C3)=O', 'C=C=CC(OCC)=O']
    # p1 = 'O=S(C1=C2C=CC=C1)(N3C2(C4=CC=CC=C4)C(C(OCC)=O)=CC3)=O'
    # p2 = 'O=S(C1=C2C=CC=C1)(N3C2(C4=CC=CC=C4)CC3=CC(OCC)=O)=O'
    # c1 = 98
    # c2 = 89
    # l1 = l2 = 314

    print('考虑反应条件')
    Predict.predict_by_smiles(True, rs, p1, c1, l1)
    Predict.predict_by_smiles(True, rs, p2, c2, l2)
    print('不考虑反应条件')
    Predict.predict_by_smiles(False, rs, p1, 0, 0)
    Predict.predict_by_smiles(False, rs, p2, 0, 0)
