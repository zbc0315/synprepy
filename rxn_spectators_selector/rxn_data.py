# -*- coding: utf-8 -*-
# @Time     : 2021/4/27 15:42
# @Author   : zbc
# @Email    : zhangbc0315@outlook.com
# @FileName : rxn_data.py
# @Software : PyCharm

import os

import json
import pandas as pd
import torch
from torch_geometric import data as gdata
from rdkit.Chem import AllChem

from nn_utils.gnn_data_utils import GnnDataUtils
from config import SSPath


class RxnData(gdata.Dataset):

    def __init__(self, role: str, data_type: str, root=None):
        self._data_type = data_type
        self._role = role
        if role == "cat":
            self._single_data_dp = SSPath.CAT_GNN_SINGLE_DATA_DP
            self._train_batch_data_dp = SSPath.TRAIN_CAT_GNN_BATCH_DATA_DP
            self._test_batch_data_dp = SSPath.TEST_CAT_GNN_BATCH_DATA_DP
            self._spec_tid_cor_fp = SSPath.CAT_TID_COR_FP
        else:
            self._single_data_dp = SSPath.SOL_GNN_SINGLE_DATA_DP
            self._train_batch_data_dp = SSPath.TRAIN_SOL_GNN_BATCH_DATA_DP
            self._test_batch_data_dp = SSPath.TEST_SOL_GNN_BATCH_DATA_DP
            self._spec_tid_cor_fp = SSPath.SOL_TID_COR_FP
        if data_type == "train":
            self._batch_data_dp = self._train_batch_data_dp
        else:
            self._batch_data_dp = self._test_batch_data_dp
        self._spec_tid_cor_df = pd.read_csv(self._spec_tid_cor_fp, sep='\t', encoding='utf-8')
        super(RxnData, self).__init__(None, None, None)

    @property
    def raw_dir(self):
        return ''

    @property
    def raw_file_names(self):
        return []

    @property
    def processed_dir(self):
        return self._single_data_dp

    @property
    def processed_file_names(self):
        return list(os.listdir(self.processed_dir))

    def download(self):
        pass

    def __len__(self):
        return self.len()

    def len(self):
        nums = {"cat": 143476,
                "sol": 634968}
        train_num = int(nums[self._role] * 0.9)
        test_num = nums[self._role] - train_num
        lens = {"train": train_num,
                "test": test_num}
        return lens[self._data_type]

    def get_fp(self, idx: int) -> (str, str):
        return os.path.join(self._single_data_dp, f"data_{self._role}_{idx}_r.pt"), \
               os.path.join(self._single_data_dp, f"data_{self._role}_{idx}_p.pt")

    def get(self, idx: int):
        r_fp, p_fp = self.get_fp(idx)
        return torch.load(r_fp), torch.load(p_fp)

    @property
    def num_node_features(self):
        return self.get(0)[0].x.size()[1]

    @property
    def num_edge_features(self):
        return self.get(0)[0].edge_attr.size()[1]

    @property
    def num_g(self):
        return self.get(0)[0].g.size().numel()

    @property
    def out_channels(self):
        return self.get(0)[0].y.size().numel()

    def get_batch_data(self):
        num = len(list(os.listdir(self._batch_data_dp)))
        for i in range(int(num/2)):
            r_batch_fp = self._get_batch_data_fp(self._batch_data_dp, i, "r")
            p_batch_fp = self._get_batch_data_fp(self._batch_data_dp, i, "p")
            yield torch.load(r_batch_fp), torch.load(p_batch_fp)

    def process(self):
        pass

    # =========================
    # -- prepare single data --

    @staticmethod
    def _get_spectator_times(smiles: str, count_df: pd.DataFrame):
        idx = count_df[count_df.smiles == smiles].index.tolist()[0]
        return count_df.loc[idx, "times"], idx

    @staticmethod
    def _get_mol_by_smileses(smileses: [str]):
        mol = None
        for i, smiles in enumerate(smileses):
            if i == 0:
                mol = AllChem.MolFromSmiles(smiles)
            else:
                mol = AllChem.CombineMols(mol, AllChem.MolFromSmiles(smiles))
        return mol

    @staticmethod
    def _get_rxn_data(reactant_mol, product_mol, y, g):
        r_data = GnnDataUtils.gnn_data_from_mol(reactant_mol)
        r_data.y = torch.tensor([y], dtype=torch.float)
        r_data.g = torch.tensor([g], dtype=torch.float)
        p_data = GnnDataUtils.gnn_data_from_mol(product_mol)
        p_data.y = torch.tensor([y], dtype=torch.float)
        p_data.g = torch.tensor([g], dtype=torch.float)
        return r_data, p_data

    @staticmethod
    def _write_all_tids(rxn_df: pd.DataFrame) -> [int]:
        tids_str = []
        tids_int = []
        for n, tid in enumerate(rxn_df["tid"]):
            if n % 50000 == 0:
                print(f"get all tids {len(tids_str)} from {n}")
            if not isinstance(tid, float):
                continue
            try:
                tid_str = str(int(tid))
            except:
                # print(tid)
                continue
            if tid_str not in tids_str:
                tids_str.append(tid_str)
                tids_int.append(int(tid_str))
        with open(SSPath.TIDS_FP, "w", encoding="utf-8")as f:
            f.write('\n'.join(tids_str))
        return tids_int

    @staticmethod
    def _read_all_tids():
        with open(SSPath.TIDS_FP, 'r', encoding="utf-8")as f:
            res = f.read()
            return [int(tid) for tid in res.split('\n')]

    @classmethod
    def _get_all_tids(cls, rxn_df: pd.DataFrame):
        if not os.path.exists(SSPath.TIDS_FP):
            cls._write_all_tids(rxn_df)
        return cls._read_all_tids()

    @classmethod
    def _get_one_hot_tid(cls, tid_float, all_tids) -> [int]:
        length = len(all_tids)
        res = [0]*(length+1)
        try:
            tid = int(tid_float)
        except:
            idx = length
            res[idx] = 1
            return res
        idx = all_tids.index(tid)
        res[idx] = 1
        return res

    @classmethod
    def _get_rxn_with_spectators(cls, role: str, spec_tid_cor_df):
        min_times = 10
        rxn_df = pd.read_csv(SSPath.RXN_WITH_SPECTATORS_FP, sep='\t', encoding="utf-8")
        spectators_count_fp = SSPath.CATS_COUNT_FP if role == "cat" else SSPath.SOLS_COUNT_FP
        spectators_count_filter_fps = {"cat": SSPath.CATS_COUNT_FILTER_FP,
                                       "sol": SSPath.SOLS_COUNT_FILTER_FP}
        spectators_count_filter_fp = spectators_count_filter_fps[role]
        spectators_df = pd.read_csv(spectators_count_fp, sep='\t', encoding="utf-8")

        # num = spectators_df["times"].count()
        filter_spectators_df = spectators_df[(spectators_df.times >= min_times) &
                                             (spectators_df.smiles != "0")]
        filter_spectators_df = filter_spectators_df.reset_index()
        filter_spectators_df.to_csv(spectators_count_filter_fp, sep='\t', index=True)
        num = filter_spectators_df["times"].count()

        # all_tids = cls._get_all_tids(rxn_df)

        for i, row in rxn_df.iterrows():
            y = [0]*num
            spectators = row[f"{role}s"]
            if not isinstance(spectators, str):
                continue
            spectators_smiles = json.loads(spectators)
            if len(spectators_smiles) > 1:
                print(f"too much spectators in {i}")
                continue
            spectator_smiles = "0" if len(spectators_smiles) == 0 else spectators_smiles[0]
            if spectator_smiles == "0":
                continue
            times, idx = cls._get_spectator_times(spectator_smiles, spectators_df)
            if times < min_times:
                continue
            times, idx = cls._get_spectator_times(spectator_smiles, filter_spectators_df)
            reactants_smiles = json.loads(row["reactants"])
            products_smiles = json.loads(row["products"])

            # tid_float = row["tid"]
            count_sids = len(eval(spec_tid_cor_df.sids[0]))
            if pd.isna(row['tid']):
                spec_tid_cor_one_hot = [0.0]*count_sids
            else:
                spec_tid_cor_queried_df = spec_tid_cor_df.query(f'tid=={int(row["tid"])}')
                spec_tid_cor_one_hot = eval(spec_tid_cor_queried_df.sids[spec_tid_cor_queried_df.index[0]])
            # tid_one_hot = cls._get_one_hot_tid(tid_float, all_tids)
            try:
                reactant_mol = cls._get_mol_by_smileses(reactants_smiles)
                product_mol = cls._get_mol_by_smileses(products_smiles)
            except Exception as e:
                print(e)
                print("="*50)
                print(i)
                continue
            if reactant_mol is None or product_mol is None:
                continue
            y[idx] = 1
            r_data, p_data = cls._get_rxn_data(reactant_mol, product_mol, y, spec_tid_cor_one_hot)
            yield row["rid"], r_data, p_data, idx, times

    def prepare_single_data(self):
        data_info = {"rid": [], "gnn_data_id": [], "spectator_id": [], "times": []}
        for n, (rid, r_data, p_data, s_idx, times) in enumerate(self._get_rxn_with_spectators(self._role, self._spec_tid_cor_df)):
            if n % 10000 == 0:
                print(n)
            data_info["rid"].append(rid)
            data_info["gnn_data_id"].append(n)
            data_info["spectator_id"].append(s_idx)
            data_info["times"].append(times)
            r_fp, p_fp = self.get_fp(n)
            torch.save(r_data, r_fp)
            torch.save(p_data, p_fp)
        data_info_df = pd.DataFrame(data_info)
        data_info_fp = SSPath.CAT_GNN_DATA_INFO_FP if self._role == "cat" else SSPath.SOL_GNN_DATA_INFO_FP
        data_info_df.to_csv(data_info_fp, sep='\t', encoding="utf-8", index=False)

    # ========================
    # -- prepare batch data --

    @classmethod
    def _get_batch_data_fp(cls, batch_dp: str, batch_id: int, r_p: str) -> str:
        return os.path.join(batch_dp, f"batch_data_{batch_id*1000}_{r_p}.pt")

    @classmethod
    def _save_data_list(cls, r_data_list: [], p_data_list: [], batch_dp: str, batch_id: int):
        r_batch_data = gdata.Batch.from_data_list(r_data_list)
        r_fp = cls._get_batch_data_fp(batch_dp, batch_id, "r")
        torch.save(r_batch_data, r_fp)

        p_batch_data = gdata.Batch.from_data_list(p_data_list)
        p_fp = cls._get_batch_data_fp(batch_dp, batch_id, "p")
        torch.save(p_batch_data, p_fp)

    @classmethod
    def prepare_batch_data(cls, role: str):
        single_data_dps = {"cat": SSPath.CAT_GNN_SINGLE_DATA_DP,
                           "sol": SSPath.SOL_GNN_SINGLE_DATA_DP}
        single_data_info_fps = {"cat": SSPath.CAT_GNN_DATA_INFO_FP,
                                "sol": SSPath.SOL_GNN_DATA_INFO_FP}
        batch_data_dps = {"cat": {"train": SSPath.TRAIN_CAT_GNN_BATCH_DATA_DP,
                                  "test": SSPath.TEST_CAT_GNN_BATCH_DATA_DP},
                          "sol": {"train": SSPath.TRAIN_SOL_GNN_BATCH_DATA_DP,
                                  "test": SSPath.TEST_SOL_GNN_BATCH_DATA_DP}}
        single_data_info_fp = single_data_info_fps[role]
        single_data_info_df = pd.read_csv(single_data_info_fp, sep='\t', encoding='utf-8')
        # single_data_info_df = single_data_info_df[single_data_info_df.times >= 10]
        single_data_info_df = single_data_info_df.sample(frac=1)
        single_data_info_df.reset_index(inplace=True)
        num = single_data_info_df["times"].count()
        train_num = int(num * 0.9)
        batch_data_dp = batch_data_dps[role]["train"]
        batch_id = 0
        r_data_list = []
        p_data_list = []
        for i, row in single_data_info_df.iterrows():
            if i == train_num:
                if len(r_data_list) != 0:
                    cls._save_data_list(r_data_list, p_data_list, batch_data_dp, batch_id)
                    r_data_list = []
                    p_data_list = []
                batch_data_dp = batch_data_dps[role]["test"]
                batch_id = 0
            if len(r_data_list) == 500:
                print(f"{batch_id} - {batch_data_dp}")
                cls._save_data_list(r_data_list, p_data_list, batch_data_dp, batch_id)
                batch_id += 1
                r_data_list = []
                p_data_list = []
            else:
                r_data = torch.load(os.path.join(single_data_dps[role], f"data_{role}_{i}_r.pt"))
                p_data = torch.load(os.path.join(single_data_dps[role], f"data_{role}_{i}_p.pt"))

                r_data_list.append(r_data)
                p_data_list.append(p_data)
        if len(r_data_list) != 0:
            cls._save_data_list(r_data_list, p_data_list, batch_data_dp, batch_id)


if __name__ == '__main__':
    # RxnData("cat", "train").prepare_single_data()
    # RxnData.prepare_batch_data("cat")
    RxnData("sol", "train").prepare_single_data()
    RxnData.prepare_batch_data("sol")
