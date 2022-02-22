#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Time    : 2022/1/10 13:29
# @Author  : zhangbc0315@outlook.com
# @File    : rct_ret_dataset.py
# @Software: PyCharm

from typing import Union, List, Tuple
import logging
import os
import random

import pandas as pd
import torch
from tqdm import tqdm
from rdkit.Chem import AllChem
import torch_geometric.data as gdata

from chem_utils.chem_handler import RxnHandler
from rxn_template_selector.filter_tids_data import FilterTidsData
from data_utils import RxnData, RxnDataType
from data_utils.rxn_template_data import RxnTemplateData
from data_utils.spectator_data import SpectatorData
from utils.gnn_data_utils import GNNDataUtils
from config.config import Config
from config.rcts_config import RCTSConfig
from config.rets_config import RETSConfig


class RCTRETDataset(gdata.InMemoryDataset):
    
    def __init__(self, train_or_test: str, rt_config: Union[RCTSConfig, RETSConfig]):
        self._config = rt_config

        if train_or_test == "train":
            self._rxn_code_and_rt_code_fp = self._config.train_rxn_code_and_rt_code_fp
            root = rt_config.train_temp_dp
        elif train_or_test == "test":
            self._rxn_code_and_rt_code_fp = self._config.test_rxn_code_and_rt_code_fp
            root = rt_config.test_temp_dp
        else:
            raise AttributeError(f"Except 'train' or 'test', but get '{train_or_test}'")

        super(RCTRETDataset, self).__init__(root)
        self.data, self.slices = torch.load(self.processed_paths[0])

    @property
    def num_global_features(self) -> int:
        return 0
        if self[0].g is None:
            return 0
        else:
            return self[0].g.size().numel()

    @property
    def raw_file_names(self) -> Union[str, List[str], Tuple]:
        return []

    @property
    def processed_file_names(self) -> Union[str, List[str], Tuple]:
        if self._config.is_retro:
            fn = f"product.{self._config.name}"
        else:
            fn = f"reactants.{self._config.name}"
        return [fn]

    def download(self):
        pass

    def process(self):
        cat_data = SpectatorData(self._config.cat_code_and_ohi_fp)
        sol_data = SpectatorData(self._config.sol_code_and_ohi_fp)
        if not os.path.exists(self._config.train_rxn_code_and_rt_code_fp):
            self._split_train_test_data()
        # rxn_code_and_template_code_df = pd.read_csv(self._rxn_code_and_rt_code_fp, sep='\t', encoding='utf-8')
        rxn_data = RxnData(self._config.rxn_data_tsv_fp, None, RxnDataType.TRAIN_TEST, self._config.evaluate_year)
        rxn_template_data = RxnTemplateData(self._config.rxn_template_tsv_fp)
        rxn_data.link_template(self._rxn_code_and_rt_code_fp)
        data_list = []
        with tqdm(total=len(rxn_data))as pbar:
            for r in rxn_data.get_all_rxn(False):
                rxn_code = r.rxn_code
                rt_code = r.rct_code if self._config.col_code == 'rct_code' else r.ret_code
                rxn_smiles = r.rxn_smiles
                if self._config.is_retro:
                    smiles = rxn_smiles.split('>')[-1]
                    g = None
                else:
                    smiles = rxn_smiles.split('>')[0]
                    g = torch.tensor([cat_data.get_one_hot_vec_by_mol_codes(r.catalysts_codes) +
                                     sol_data.get_one_hot_vec_by_mol_codes(r.solvents_codes)], dtype=torch.float)
                data = GNNDataUtils.get_gnn_data_from_smiles(smiles)
                data.g = g
                data.y = torch.tensor([rxn_template_data.get_one_hot_idx_by_rt_code(rt_code)], dtype=torch.long)
                data_list.append(data)

                pbar.update(1)
        data, slices = self.collate(data_list)
        torch.save((data, slices), self.processed_paths[0])

    def _split_train_test_data(self):
        logging.info(f"Splitting train and test to {self._config.train_rxn_code_and_rt_code_fp} and {self._config.test_rxn_code_and_rt_code_fp}")

        rxn_code_and_rt_code_df = pd.read_csv(self._config.rxn_code_rt_code_tsv_fp, sep='\t', encoding='utf-8')
        filter_rxn_code_and_rt_code_df = rxn_code_and_rt_code_df.query(f"{self._config.col_num}>={self._config.min_num_covered_rxns_by_rxn_template}")
        filter_rxn_code_and_rt_code_df = filter_rxn_code_and_rt_code_df.sample(frac=1)
        num_train = int(len(filter_rxn_code_and_rt_code_df) * 0.9)
        filter_rxn_code_and_rt_code_df.iloc[:num_train].to_csv(self._config.train_rxn_code_and_rt_code_fp, sep='\t', encoding='utf-8', index=False)
        filter_rxn_code_and_rt_code_df.iloc[num_train:].to_csv(self._config.test_rxn_code_and_rt_code_fp, sep='\t', encoding='utf-8', index=False)


if __name__ == "__main__":
    config = Config("../config.json")
    RCTRETDataset("train", config.rets_config)
    RCTRETDataset("test", config.rets_config)
