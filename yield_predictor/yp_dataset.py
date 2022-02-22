#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Time    : 2022/2/12 15:20
# @Author  : zhangbc0315@outlook.com
# @File    : yp_dataset.py
# @Software: PyCharm

from typing import Union, List, Tuple
import os
import random

import torch
from tqdm import tqdm
import pandas as pd
import torch_geometric.data as gdata

from utils import GNNDataUtils
from config import YPConfig, Config
from data_utils import Mol, MolData, Rxn, RxnData, RxnDataType, SpectatorData, SpectatorType
from yield_predictor.data_analyzer import DataAnalyzer


class YPDataset(gdata.InMemoryDataset):

    def __init__(self, train_or_test: str, yp_config: YPConfig, r_or_p: str):
        self._train_or_test = train_or_test
        self._yp_config = yp_config
        self._rxn_data, self._rxn_codes = self._get_rxn_df()
        self._cat_data = SpectatorData(yp_config.cat_code_to_ohi_fp)
        self._sol_data = SpectatorData(yp_config.sol_code_to_ohi_fp)
        self._r_or_p = r_or_p

        if train_or_test == "train":
            root = yp_config.train_temp_dp
        elif train_or_test == "test":
            root = yp_config.test_temp_dp
        else:
            raise AttributeError(f"Except 'train' or 'test', but get '{train_or_test}'")
        super(YPDataset, self).__init__(root)
        if r_or_p == 'r':
            self.data, self.slices = torch.load(self.processed_paths[0])
        else:
            self.data, self.slices = torch.load(self.processed_paths[1])

    @property
    def num_global_features(self) -> int:
        return self[0].g.size().numel()

    @property
    def raw_file_names(self) -> Union[str, List[str], Tuple]:
        return []

    @property
    def processed_file_names(self) -> List[str]:
        return ['rxn.yield.reactants', 'rxn.yield.products']

    def download(self):
        pass

    def process(self):
        data_list = []
        with tqdm(total=len(self._rxn_codes))as pbar:
            pbar.set_description(f"{self._train_or_test} {self._r_or_p}")
            for i, rxn_code in enumerate(self._rxn_codes):
            # for i, row in self._rxn_df.iterrows():
            #     if row.rxn_code not in self._rxn_codes:
            #         continue
                pbar.update(1)
                rxn_d = self._rxn_data.get_rxn_by_rxn_code(rxn_code, False)
                rxn_smiles = rxn_d.rxn_smiles
                if self._r_or_p == 'r':
                    smiles = rxn_smiles.split('>')[0]
                    data = GNNDataUtils.get_gnn_data_from_smiles(smiles)
                    cat_ohv = self._cat_data.get_one_hot_vec_by_mol_codes(rxn_d.catalysts_codes)
                    sol_ohv = self._sol_data.get_one_hot_vec_by_mol_codes(rxn_d.solvents_codes)
                    data.g = torch.tensor([cat_ohv+sol_ohv], dtype=torch.float)
                else:
                    smiles = rxn_smiles.split('>')[-1]
                    data = GNNDataUtils.get_gnn_data_from_smiles(smiles)
                    data.y = torch.tensor([rxn_d.yield_product], dtype=torch.float)
                data.key = rxn_d.rxn_code
                # cat_ohv = self._cat_data.get_one_hot_vec_by_mol_codes(eval(row.catalysts_codes))
                # sol_ohv = self._sol_data.get_one_hot_vec_by_mol_codes(eval(row.solvents_codes))
                # yield_product = row.yield_product
                # data = GNNDataUtils.get_gnn_data_from_rxn_smiles(rxn_smiles)
                # data.g = torch.tensor([cat_ohv+sol_ohv], dtype=torch.float)
                # data.y = yield_product
                data_list.append(data)
        data, slices = self.collate(data_list)
        if self._r_or_p == 'r':
            torch.save((data, slices), self.processed_paths[0])
        else:
            torch.save((data, slices), self.processed_paths[1])

    def _split_train_test(self, rxn_codes: [str]):
        random.shuffle(rxn_codes)
        train_num = int(len(rxn_codes) * 0.9)
        train_codes = rxn_codes[:train_num]
        test_codes = rxn_codes[train_num:]
        with open(self._yp_config.train_rxn_codes_fp, 'w', encoding='utf-8')as f:
            f.write('\n'.join(train_codes))
        with open(self._yp_config.test_rxn_codes_fp, 'w', encoding='utf-8')as f:
            f.write('\n'.join(test_codes))

    @classmethod
    def _get_rxn_codes(cls, rxn_codes_fp: str) -> [str]:
        rxn_codes = []
        with open(rxn_codes_fp, 'r', encoding='utf-8')as f:
            for line in f.readlines():
                line = line.strip()
                if len(line) == 0:
                    continue
                rxn_codes.append(line)
        return rxn_codes

    def _get_rxn_df(self) -> (pd.DataFrame, [str]):
        rxn_data = RxnData(self._yp_config.rxn_data_fp,
                           None,
                           RxnDataType.TRAIN_TEST,
                           self._yp_config.evaluate_year)
        rxn_data.init_yield()
        if not os.path.exists(self._yp_config.train_rxn_codes_fp) or not os.path.exists(self._yp_config.test_rxn_codes_fp):
            rxn_codes = list(rxn_data.rxn_df.rxn_code)
            # rxn_codes = DataAnalyzer.duplicate_rare_data(rxn_data.rxn_df)
            self._split_train_test(rxn_codes)
        if self._train_or_test == 'train':
            rxn_codes = self._get_rxn_codes(self._yp_config.train_rxn_codes_fp)
        elif self._train_or_test == 'test':
            rxn_codes = self._get_rxn_codes(self._yp_config.test_rxn_codes_fp)
        else:
            raise ValueError(f"expect 'train' or 'test', but get {self._train_or_test}")
        return rxn_data, rxn_codes


if __name__ == "__main__":
    pass
    # config = Config("../config.json")
    # YPDataset("test", config.yp_config)
    # YPDataset("train", config.yp_config)
