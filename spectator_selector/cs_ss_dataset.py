#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Time    : 2022/2/21 13:40
# @Author  : zhangbc0315@outlook.com
# @File    : cs_ss_dataset.py
# @Software: PyCharm
from typing import Union, List, Tuple
import os

import pandas as pd
import torch
from torch_geometric import data as gdata
from tqdm import tqdm

from config.cs_ss_config import CSSSConfig
from data_utils.spectator_data import SpectatorData
from data_utils import RxnData, RxnDataType
from utils import GNNDataUtils, ModelUtils


class SSDataset(gdata.InMemoryDataset):

    def __init__(self, train_or_test: str, config: CSSSConfig, r_or_p: str):
        self._config = config
        self._r_or_p = r_or_p

        if train_or_test == "train":
            self._train_or_test_data_fp = self._config.train_rxn_code_fp
            root = self._config.train_temp_dp
        elif train_or_test == "test":
            self._train_or_test_data_fp = self._config.test_rxn_code_fp
            root = self._config.test_temp_dp
        else:
            raise AttributeError(f"Except 'train' or 'test', but get '{train_or_test}'")
        super(SSDataset, self).__init__(root)
        if r_or_p == 'r':
            self.data, self.slices = torch.load(self.processed_paths[0])
        else:
            self.data, self.slices = torch.load(self.processed_paths[1])

    @property
    def processed_file_names(self) -> List[str]:
        return ['rxn.spectator.reactants', 'rxn.epectator.products']

    @property
    def raw_file_names(self) -> Union[str, List[str], Tuple]:
        return []

    def download(self):
        pass

    def process(self):
        spectator_data = SpectatorData(self._config.spectator_code_to_ohi_fp)
        if not os.path.exists(self._config.train_temp_dp):
            self._split_data()
        train_or_test_data_df = pd.read_csv(self._train_or_test_data_fp, sep='\t', encoding='utf-8')
        data_list = []
        with tqdm(total=len(train_or_test_data_df))as pbar:
            pbar.set_description(f"process spectator for {self._r_or_p}")
            for i, row in train_or_test_data_df.iterrows():
                rxn_smiles = row.rxn_smiles
                if self._r_or_p == 'r':
                    smiles = rxn_smiles.split('>')[0]
                else:
                    smiles = rxn_smiles.split('>')[-1]
                data = GNNDataUtils.get_gnn_data_from_smiles(smiles)
                data.y = torch.tensor([spectator_data.get_one_hot_idx_by_mol_code(row.spectator_code)], dtype=torch.long)
                data_list.append(data)
                pbar.update(1)
        data, slices = self.collate(data_list)
        if self._r_or_p == 'r':
            torch.save((data, slices), self.processed_paths[0])
        else:
            torch.save((data, slices), self.processed_paths[1])

    def _split_data(self):
        rxn_code_and_spe_code_df = pd.read_csv(self._config.rxn_code_and_spe_code_fp)
        rxn_code_and_spe_code_df = rxn_code_and_spe_code_df.sample(frac=1)
        num_train = int(len(rxn_code_and_spe_code_df)*0.9)
        rxn_code_and_spe_code_df.iloc[:num_train].to_csv(self._config.train_rxn_code_fp, sep='\t', encoding='utf-8', index=False)
        rxn_code_and_spe_code_df.iloc[num_train:].to_csv(self._config.test_rxn_code_fp, sep='\t', encoding='utf-8', index=False)


if __name__ == "__main__":
    pass
