#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Time    : 2022/2/12 10:34
# @Author  : zhangbc0315@outlook.com
# @File    : spectators_prepare.py
# @Software: PyCharm

from tqdm import tqdm
import pandas as pd

from config.config import Config
from data_utils.rxn_data import MolData, Rxn, RxnData, RxnDataType


class SpectatorsPrepare:

    def __init__(self):
        self.config = Config("config.json")
        # self.mol_data = MolData(self.config.mol_data_tsv_file_path)
        self.rxn_data = RxnData(self.config.rxn_data_tsv_file_path,
                                None,
                                RxnDataType.TRAIN_TEST,
                                self.config.evaluate_year)

        self.cat_to_count = {}
        self.sol_to_count = {}

        self.num_cat_and_sol = 0  # 具有一个常见催化剂和一个常见溶剂的反应数目
        self.num_cat_and_none_sol = 0
        self.num_cat_and_other_sol = 0
        self.num_cat_and_sols = 0

        self.num_none_cat_and_sol = 0
        self.num_none_cat_and_none_sol = 0
        self.num_none_cat_and_other_sol = 0
        self.num_none_cat_and_sols = 0

        self.num_other_cat_and_sol = 0
        self.num_other_cat_and_none_sol = 0
        self.num_other_cat_and_other_sol = 0
        self.num_other_cat_and_sols = 0

        self.num_cats_and_sol = 0
        self.num_cats_and_none_sol = 0
        self.num_cats_and_other_sol = 0
        self.num_cats_and_sols = 0

    def count_spectators(self):
        with tqdm(total=len(self.rxn_data))as pbar:
            for n, rxn in enumerate(self.rxn_data.get_all_rxn(need_mol_inchi_smiles=False)):
                pbar.set_description("count spectators")
                pbar.update(1)
                cats_codes = rxn.catalysts_codes
                sols_codes = rxn.solvents_codes
                if len(cats_codes) == 1:
                    cat_code = cats_codes[0]
                    if cat_code in self.cat_to_count.keys():
                        self.cat_to_count[cat_code] += 1
                    else:
                        self.cat_to_count[cat_code] = 1
                if len(sols_codes) == 1:
                    sol_code = sols_codes[0]
                    if sol_code in self.sol_to_count.keys():
                        self.sol_to_count[sol_code] += 1
                    else:
                        self.sol_to_count[sol_code] = 1

    def dump_spectators(self, spectator_to_count: {str, int}, tsv_fp: str, min_count: int,
                        filter_rxn_fp: str, cat_or_sol: str):
        mol_code_and_one_hot_idx = {"mol_code": [], "one_hot_idx": []}
        one_hot_idx = 0
        mol_codes = set()
        for mol_code, count in spectator_to_count.items():
            if count >= min_count:
                mol_code_and_one_hot_idx["mol_code"].append(mol_code)
                mol_code_and_one_hot_idx["one_hot_idx"].append(one_hot_idx)
                one_hot_idx += 1
                mol_codes.add(mol_code)
        df = pd.DataFrame(mol_code_and_one_hot_idx)
        df.to_csv(tsv_fp, sep='\t', encoding='utf-8', index=False)

        rxn_df = pd.DataFrame(columns=["rxn_code", "rxn_smiles", "spectator_code"])
        with tqdm(total=len(self.rxn_data)) as pbar:
            for n, rxn in enumerate(self.rxn_data.get_all_rxn(False)):
                pbar.set_description("dump spectators")
                pbar.update(1)
                if cat_or_sol == 'cat':
                    spec_codes = rxn.catalysts_codes
                elif cat_or_sol == 'sol':
                    spec_codes = rxn.solvents_codes
                else:
                    raise AttributeError(f"cat_or_sol should be 'cat' or 'sol', but get {cat_or_sol}")
                if len(spec_codes) != 1 or spec_codes[0] not in mol_codes:
                    continue
                rxn_df = rxn_df.append({"rxn_code": rxn.rxn_code, "rxn_smiles": rxn.rxn_smiles, "spectator_code": spec_codes[0]}, ignore_index=True)
        rxn_df.to_csv(filter_rxn_fp, sep='\t', encoding='utf-8', index=False)

    def process(self):
        self.count_spectators()
        self.dump_spectators(self.cat_to_count,
                             self.config.cat_code_and_one_hot_idx_tsv_file_path,
                             self.config.min_num_covered_rxns_by_cat,
                             self.config.rxn_code_and_cat_code_fp,
                             "cat")
        self.dump_spectators(self.sol_to_count,
                             self.config.sol_code_and_one_hot_idx_tsv_file_path,
                             self.config.min_num_covered_rxns_by_sol,
                             self.config.rxn_code_and_sol_code_fp,
                             "sol")


if __name__ == "__main__":
    SpectatorsPrepare().process()
