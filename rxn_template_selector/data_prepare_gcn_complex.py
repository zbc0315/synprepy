#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Time    : 2021/8/26 14:23
# @Author  : zhangbc0315@outlook.com
# @File    : data_prepare_gcn_complex.py
# @Software: PyCharm
import os
from random import shuffle

import torch
from torch_geometric.data import Batch

from data_utils.tid_utils import TidUtils
from nn_utils.gnn_data_utils import GnnDataUtils
from config import CTPath
from rmp_database.old_chem_db import OldChemDb


class DataPrepareGnnComplex:

    _odb = OldChemDb()

    @classmethod
    def get_tid_data(cls):
        odb = OldChemDb()
        for tid_data in odb.get_data_iter('public.rxn_template_10_filter', ['tid', 'rids'], None):
            rids = [int(rid) for rid in tid_data['rids'].split('.')]
            yield tid_data['tid'], rids

    @classmethod
    def get_p_inchi(cls, rid):
        rxn_data = cls._odb.search_one_data('public.cleaned_rxn', ['rxn_code'], f'rid={rid}')
        if rxn_data is None:
            return None

        p_ids = [int(pid) for pid in rxn_data['rxn_code'].split('>')[-1].split('.')]
        if len(p_ids) != 1:
            return None

        mol_data = cls._odb.search_one_data('public.mols', ['inchi'], f'mid={p_ids[0]}')
        if mol_data is None:
            return None
        return mol_data['inchi']

    @classmethod
    def get_all_rxn_data(cls):
        for tid, rids in cls.get_tid_data():
            for rid in rids:
                p_inchi = cls.get_p_inchi(rid)
                yield p_inchi, tid

    @classmethod
    def prepare_data(cls):
        tid_to_idx = TidUtils.load_comp_tid_to_idx()
        num = 0
        for p_inchi, tid in cls.get_all_rxn_data():
            if p_inchi is None:
                print('can not get product`s inchi')
                continue
            data = GnnDataUtils.gnn_data_from_inchi(p_inchi)
            if data is None:
                print(f'can not get gnn data from {p_inchi}')
                continue

            y = [0]*len(tid_to_idx.keys())
            y[tid_to_idx[tid]] = 1
            data.y = torch.tensor([y], dtype=torch.float)
            torch.save(data, os.path.join(CTPath.GCN_SINGLE_DATA_DP, f'{num}.pt'))
            num += 1

    # =========================

    @classmethod
    def get_all_tids(cls):
        tids = set()
        odb = OldChemDb()
        for tid_data in odb.get_data_iter('public.rxn_template_10_filter', ['tid'], None):
            tids.add(tid_data['tid'])

        tids = list(tids)
        tids = sorted(tids, key=lambda x: x)
        return tids

    @classmethod
    def prepare_tids(cls):
        tids = cls.get_all_tids()
        TidUtils.save_comp_tids(tids)

    # ===========================

    @classmethod
    def process(cls):
        if not os.path.exists(CTPath.COMP_TIDS_FP):
            cls.prepare_tids()
        cls.prepare_data()

    # ===========================

    @classmethod
    def get_train_test_fns(cls):
        fns = list(os.listdir(CTPath.GCN_SINGLE_DATA_DP))
        length = len(fns)
        shuffle(fns)
        return fns[:int(0.9*length)], fns[int(0.9*length):]

    @classmethod
    def save_batch_data(cls, data_list: [], idx: int, dp: str):
        fp = os.path.join(dp, f'{idx*1000}-{idx*1000 + len(data_list) - 1}.pt')
        print(fp)
        batch = Batch.from_data_list(data_list)
        torch.save(batch, fp)

    @classmethod
    def get_batch_data(cls, fns: [str], dp: str):
        data_list = []
        num = 0
        for fn in fns:
            if len(data_list) == 1000:
                cls.save_batch_data(data_list, num, dp)
                num += 1
                data_list = []
            data = torch.load(os.path.join(CTPath.GCN_SINGLE_DATA_DP, fn))
            data_list.append(data)
        if len(data_list) != 0:
            cls.save_batch_data(data_list, num, dp)

    @classmethod
    def process_batch(cls):
        train_fns, test_fns = cls.get_train_test_fns()
        cls.get_batch_data(train_fns, CTPath.GCN_TRAIN_BATCH_DATA_DP)
        cls.get_batch_data(test_fns, CTPath.GCN_TEST_BATCH_DATA_DP)


if __name__ == "__main__":
    DataPrepareGnnComplex.process()
    DataPrepareGnnComplex.process_batch()
