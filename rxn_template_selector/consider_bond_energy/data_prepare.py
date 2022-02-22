#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Time    : 2021/8/18 16:44
# @Author  : zhangbc0315@outlook.com
# @File    : data_prepare.py
# @Software: PyCharm
import os
import random
from datetime import datetime
import sys
import time
from multiprocessing import Pool

import torch
from torch_geometric.data import Batch

from config import BTPath
from data_utils.data_utils import DataUtils
from bond_energy_predictor.bep_data_utils import BepDataUtils


class DataPrepare:

    data_idx = 0
    start = datetime.now()

    # single_data_dps = {"train": BTPath.CBE_SINGLE_DATA_TRAIN_DP,
    #                    "test": BTPath.CBE_SINGLE_DATA_TEST_DP}

    batch_data_dps = {"train": BTPath.CBE_BATCH_DATA_TRAIN_DP,
                      "test": BTPath.CBE_BATCH_DATA_TEST_DP}

    c_path = 'C:\\Users\\zhang\\Documents\\Data\\RxnMolPredictor\\BaseTemp\\consider_bond_energy'
    single_data_dps = {"train": os.path.join(c_path, 'train_single_data'),
                       "test": os.path.join(c_path, 'test_single_data')}
    # batch_data_dps = {'train': os.path.join(c_path, 'train_batch_1000_data'),
    #                   'test': os.path.join(c_path, 'test_batch_1000_data')}

    @classmethod
    def load_saved_idx(cls, train_test: str):
        cls.data_idx = len(list(os.listdir(cls.single_data_dps[train_test])))

    @classmethod
    def prepare_many_data(cls, idxes_smileses_ys, single_dp):
        smileses = []
        for (_, smiles, __) in idxes_smileses_ys:
            smileses.append(smiles)
        ys
        data_list = BepDataUtils.get_data_list_by_smileses(smileses)
        for data in data_list:
            data_fp = os.path.join(single_dp, f"{cls.data_idx}.pt")
            torch.save(data, data_fp)
            cls.data_idx += 1
            if cls.data_idx % 1000 == 0:
                end = datetime.now()
                print(f"{cls.data_idx}: {(end - cls.start).seconds}")
                cls.start = datetime.now()

    @classmethod
    def prepare_train_test(cls, train_test: str, pid: int):
        # data_tn = f"ml_data.{train_test}_gnn"
        cls.data_idx = 0
        cls.load_saved_idx(train_test)
        data_tn = f"C:/Users/zhang/OneDrive/Projects/RXN-PREDICTOR/Data/DataBase/rxndb_ml_data_{train_test}_gnn.tsv"
        single_dp = cls.single_data_dps[train_test]

        n = 0
        idxes_smileses_ys = []
        with torch.no_grad():
            for idx, smiles, encoded_tid in DataUtils.iter_db_data(data_tn, pid, cls.data_idx):
                idxes_smileses_ys.append((idx, smiles, encoded_tid))
                if len(idxes_smileses_ys) == 150:
                    cls.prepare_many_data(idxes_smileses_ys, single_dp)
                    idxes_smileses_ys = []
            if len(idxes_smileses_ys) > 0:
                cls.prepare_many_data(idxes_smileses_ys, single_dp)

        # with torch.no_grad():
        #     for idx, smiles, encoded_tid in DataUtils.iter_db_data(data_tn, pid):
        #         fp = os.path.join(single_dp, f"{idx}.pt")
        #         try:
        #             data = BepDataUtils.get_data_by_smiles(smiles)
        #             data.key = smiles
        #         except Exception as e:
        #             print(smiles)
        #             raise e
        #         data.y = torch.tensor([encoded_tid], dtype=torch.float)
        #         torch.save(data, fp)
        #         n += 1
        #         if n % 1000 == 0:
        #             end = datetime.now()
        #             print(f"{n}: {(end - cls.start).seconds}")
        #             cls.start = datetime.now()

    @classmethod
    def process(cls, pid):
        time.sleep(pid*5)
        print(f"start process: {pid}")
        cls.prepare_train_test("train", pid)
        cls.prepare_train_test("test", pid)

    @classmethod
    def save_batch(cls, batch_id, data_list, batch_dp):
        print(f"{batch_id} - {len(data_list)}")
        batch_fn = f"{batch_id*1000}-{batch_id*1000 + len(data_list) - 1}.pt"
        batch_fp = os.path.join(batch_dp, batch_fn)
        batch = Batch.from_data_list(data_list)
        torch.save(batch, batch_fp)

    @classmethod
    def prepare_train_test_batch(cls, train_test: str):
        single_dp = cls.single_data_dps[train_test]
        batch_dp = cls.batch_data_dps[train_test]
        fns = os.listdir(single_dp)
        random.shuffle(fns)

        data_list = []
        batch_id = 0
        for i, fn in enumerate(fns):
            if len(data_list) == 1000:
                cls.save_batch(batch_id, data_list, batch_dp)
                batch_id += 1
                data_list = []
            data = torch.load(os.path.join(single_dp, fn))
            if len(data.edge_attr.shape) < 2:
                continue
            elif data.edge_attr.shape[1] == 13:
                continue
            data_list.append(torch.load(os.path.join(single_dp, fn)))
        if len(data_list) > 0:
            cls.save_batch(batch_id, data_list, batch_dp)

    @classmethod
    def process_batch(cls):
        cls.prepare_train_test_batch('train')
        cls.prepare_train_test_batch('test')


if __name__ == "__main__":
    # num_ps = 10
    # p = Pool(num_ps)
    # for i in range(num_ps):
    #     p.apply_async(DataPrepare.process, args=(i, ))
    # p.close()
    # p.join()
    # DataPrepare.process(0)
    DataPrepare.process_batch()
