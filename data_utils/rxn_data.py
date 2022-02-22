#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Time    : 2022/2/11 14:33
# @Author  : zhangbc0315@outlook.com
# @File    : rxn_data.py
# @Software: PyCharm

from typing import Optional, Iterator
from enum import Enum

import pandas as pd

from data_utils.rxn import Rxn
from data_utils.mol_data import MolData


class RxnDataType(Enum):

    TRAIN_TEST = 0,
    EVAL = 1


class RxnData:

    def __init__(self, rxn_fp: str, mol_data: Optional[MolData], data_type: RxnDataType, eval_year: int):
        rxn_df = pd.read_csv(rxn_fp, sep='\t', encoding='utf-8')
        self._rxn_df = self._get_rxn_by_data_type(rxn_df, data_type, eval_year)
        self._mol_data = mol_data

    def __len__(self):
        return len(self._rxn_df)

    @classmethod
    def _get_rxn_by_data_type(cls, rxn_df: pd.DataFrame, data_type: RxnDataType, eval_year: int):
        if data_type == RxnDataType.TRAIN_TEST:
            return rxn_df.query(f"create_year!={eval_year}")
        elif data_type == RxnDataType.EVAL:
            return rxn_df.query(f"create_year=={eval_year}")
        else:
            raise ValueError(f"需要DataType格式的data_type，但得到的是格式为{type(data_type)}的{data_type}")

    def init_yield(self, min_yield: int = 0, max_yield: int = 100):
        self._rxn_df = self._rxn_df.dropna(subset=["yield_product"])
        self._rxn_df = self._rxn_df.query(f"yield_product >= {min_yield}")
        self._rxn_df = self._rxn_df.query(f"yield_product <= {max_yield}")

    def link_template(self, rxn_code_and_rt_code_fp: str):
        if 'rct_code' not in self._rxn_df.columns:
            rxn_code_rt_code_df = pd.read_csv(rxn_code_and_rt_code_fp, sep='\t', encoding='utf-8')
            # self._rxn_df = pd.concat([self._rxn_df, rxn_code_rt_code_df], axis=1, join='inner')
            self._rxn_df = pd.merge(self._rxn_df, rxn_code_rt_code_df, on=['rxn_code'])

    def init_template_filter(self, col_name: str, min_num: int):
        self._rxn_df = self._rxn_df.quey(f"{col_name}>={min_num}")

    @property
    def rxn_df(self) -> pd.DataFrame:
        return self._rxn_df

    def get_all_rxn(self, need_mol_inchi_smiles: bool) -> Iterator[Rxn]:
        for _, row in self._rxn_df.iterrows():
            reactants_codes = eval(row["reactants_codes"])
            catalysts_codes = eval(row["catalysts_codes"])
            solvents_codes = eval(row["solvents_codes"])
            rxn = Rxn(row["rxn_code"],
                      row["rxn_mols_code"],
                      reactants_codes,
                      row["product_code"],
                      catalysts_codes,
                      solvents_codes,
                      row["yield_product"],
                      row["rxn_smiles"],
                      row["create_year"])
            if "rct_code" in self._rxn_df.columns:
                rxn.rct_code = row["rct_code"]
                rxn.ret_code = row["ret_code"]
            if need_mol_inchi_smiles:
                rxn.init_rxn(self._mol_data)
            yield rxn

    def get_rxn_by_rxn_code(self, rxn_code: str, need_mol_inchi_smiles: bool) -> Optional[Rxn]:
        rxn_df = self._rxn_df.query(f"rxn_code=='{rxn_code}'")
        if len(rxn_df) == 0:
            return None
        rxn_df = rxn_df.reset_index()
        reactants_codes = eval(rxn_df["reactants_codes"][0])
        catalysts_codes = eval(rxn_df["catalysts_codes"][0])
        solvents_codes = eval(rxn_df["solvents_codes"][0])
        rxn = Rxn(rxn_df["rxn_code"][0],
                  rxn_df["rxn_mols_code"][0],
                  reactants_codes,
                  rxn_df["product_code"][0],
                  catalysts_codes,
                  solvents_codes,
                  rxn_df["yield_product"][0],
                  rxn_df["rxn_smiles"][0],
                  rxn_df["create_year"][0])
        if need_mol_inchi_smiles:
            rxn.init_rxn(self._mol_data)
        return rxn


if __name__ == "__main__":
    from config import Config
    from data_utils.rxn_template_data import RxnTemplateData
    config = Config('../config.json')
    rd = RxnData(config.rxn_data_tsv_file_path, None, RxnDataType.TRAIN_TEST, config.evaluate_year)
    rd.link_template(config.rxn_code_with_rt_code_tsv_file_path)
    td = RxnTemplateData(config.rxn_centralized_template_tsv_file_path)
    for r in rd.get_all_rxn(False):
        r_smi = r.rxn_smiles
        t_smi = td.get_rt_smiles_by_rt_code(r.rct_code)
        print(r_smi)
        print(t_smi)
        print(r.rct_code)
        print(td.get_one_hot_idx_by_rt_code(r.rct_code))
        print('#'*100)
