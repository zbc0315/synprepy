#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Time    : 2021/8/2 15:22
# @Author  : zhangbc0315@outlook.com
# @File    : bep_data_utils.py
# @Software: PyCharm

from typing import Tuple, List
import math

import torch
from rdkit.Chem import AllChem
import torch_geometric.data as gdata
from torch_geometric.data import Data, Batch

from bond_energy_predictor.model_gcn import Net
from nn_utils.gnn_data_utils import GnnDataUtils
from config import BTPath


class BepDataUtils:

    _symbol_to_idx = {'c': 0,
                      's': 1,
                      'o': 2,
                      'n': 3,
                      'cl': 4,
                      'h': 5,
                      'f': 6,
                      'i': 7,
                      'p': 8,
                      'br': 9,
                      'si': 10}
    _device = 'cuda'
    _bond_energy_model = Net.get_model(_device)

    @classmethod
    def get_bond_length(cls, bond, positions):
        begin_aid = bond.GetBeginAtomIdx()
        end_aid = bond.GetEndAtomIdx()
        p1 = positions[begin_aid]
        p2 = positions[end_aid]
        l2 = (p2 - p1) ** 2
        return math.sqrt(l2.sum())

    @classmethod
    def atom_to_feature(cls, atom):
        symbol = atom.GetSymbol().lower()
        res = [0]*11
        idx = cls._symbol_to_idx[symbol]
        res[idx] = 1
        return res

    @classmethod
    def get_nodes_feature(cls, rdmol):
        res = []
        for atom in rdmol.GetAtoms():
            res.append(cls.atom_to_feature(atom))
        return res

    @classmethod
    def get_new_bond_aids(cls, bond, old_aid_to_new_aid):
        res1 = [old_aid_to_new_aid[bond.GetBeginAtomIdx()], old_aid_to_new_aid[bond.GetEndAtomIdx()]]
        res2 = [old_aid_to_new_aid[bond.GetEndAtomIdx()], old_aid_to_new_aid[bond.GetBeginAtomIdx()]]
        return res1, res2

    @classmethod
    def get_connect_and_features(cls, rdmol, positions) -> Tuple[List[List[int]], List[List[float]]]:
        connect = []
        edge_features = []
        for bond in rdmol.GetBonds():
            bond_length = cls.get_bond_length(bond, positions)

            connect.append([bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()])
            connect.append([bond.GetEndAtomIdx(), bond.GetBeginAtomIdx()])

            edge_features.append([1/bond_length, 1/bond_length])
            edge_features.append([1/bond_length, 1/bond_length])

        return connect, edge_features

    @classmethod
    def get_bond_types(cls, rdmol, concern_bond) -> List[int]:
        edge_types = []
        for bond in rdmol.GetBonds():
            bond_type = 0 if bond.GetIdx() == concern_bond.GetIdx() else 1
            edge_types.append(bond_type)
            edge_types.append(bond_type)
        return edge_types

    @classmethod
    def calc_coordinates(cls, rdmol):
        rdmol = AllChem.AddHs(rdmol)
        status = AllChem.EmbedMolecule(rdmol)
        status = AllChem.UFFOptimizeMolecule(rdmol)
        conformer = rdmol.GetConformer()
        positions = conformer.GetPositions()
        return rdmol, positions

    # ====================================

    @classmethod
    def add_bond_energy_to_features(cls, bond_energies, bond_features, old_bid_to_new_bid):
        for i, be in enumerate(bond_energies):
            if i not in old_bid_to_new_bid.keys():
                continue
            new_bid = old_bid_to_new_bid[i]
            # bond_features[new_bid*2].append(be/100)
            # bond_features[new_bid*2+1].append(be/100)
            bond_features[new_bid * 2][0] = be / 100
            bond_features[new_bid * 2 + 1][0] = be / 100
        return bond_features

    @classmethod
    def get_rts_data_by_rdmol_and_100(cls, rdmol):
        x = GnnDataUtils.nodes_features_from_mol(rdmol, '')
        connect = GnnDataUtils.connect_info_from_mol(rdmol)
        edge_feature = GnnDataUtils.edges_features_from_mol(rdmol, '')
        [e.append(1) for e in edge_feature]
        data = Data(x=torch.tensor(x, dtype=torch.float),
                    edge_index=torch.tensor(connect, dtype=torch.long).t().contiguous(),
                    edge_attr=torch.tensor(edge_feature, dtype=torch.float))
        return data

    @classmethod
    def get_rts_data_by_rdmol_and_bes(cls, rdmol, bond_energies, predict_type=''):
        x, old_aid_to_new_aid = GnnDataUtils.nodes_features_from_mol(rdmol, predict_type, no_h=True, transform_ex_h_to_im_h=True)
        connect = GnnDataUtils.connect_info_from_mol(rdmol, no_h=True, old_aid_to_new_aid=old_aid_to_new_aid)
        edge_features, old_bid_to_new_bid = GnnDataUtils.edges_features_from_mol(rdmol, predict_type, no_h=True)
        edge_features = cls.add_bond_energy_to_features(bond_energies, edge_features, old_bid_to_new_bid)
        data = Data(x=torch.tensor(x, dtype=torch.float),
                    edge_index=torch.tensor(connect, dtype=torch.long).t().contiguous(),
                    edge_attr=torch.tensor(edge_features, dtype=torch.float))
        return data

    # =======================================

    @classmethod
    def get_rts_data_by_rdmol(cls, rdmol):
        rdmol, positions = cls.calc_coordinates(rdmol)
        data_list = []
        x = cls.get_nodes_feature(rdmol)
        edge_index, edge_attr = cls.get_connect_and_features(rdmol, positions)

        for bond in rdmol.GetBonds():
            edge_type = cls.get_bond_types(rdmol, bond)
            data = Data(x=torch.tensor(x, dtype=torch.float),
                        edge_index=torch.tensor(edge_index, dtype=torch.long).t().contiguous(),
                        edge_attr=torch.tensor(edge_attr, dtype=torch.float),
                        edge_type=torch.tensor(edge_type, dtype=torch.long),
                        g=torch.tensor([[1]], dtype=torch.float))
            data_list.append(data)
        batch_data = Batch.from_data_list(data_list).to(cls._device)
        bond_energies = cls._bond_energy_model(batch_data, len(data_list))

        rts_data = cls.get_rts_data_by_rdmol_and_bes(rdmol, bond_energies.view(1, len(data_list))[0].tolist())
        return rts_data

    @classmethod
    def get_data_by_smiles(cls, smiles: str):

        rdmol = AllChem.MolFromSmiles(smiles)
        try:
            data = cls.get_rts_data_by_rdmol(rdmol)
        except:
            data = cls.get_rts_data_by_rdmol_and_100(rdmol)
        # data1 = cls.get_rts_data_by_rdmol(rdmol)
        # batch = Batch.from_data_list([data, data1])
        # ea = batch.edge_attr[:, :-1]
        return data

    @classmethod
    def get_data_list_to_bde(cls, rdmol):
        rdmol, positions = cls.calc_coordinates(rdmol)
        data_list = []
        x = cls.get_nodes_feature(rdmol)
        edge_index, edge_attr = cls.get_connect_and_features(rdmol, positions)

        for bond in rdmol.GetBonds():
            edge_type = cls.get_bond_types(rdmol, bond)
            data = Data(x=torch.tensor(x, dtype=torch.float),
                        edge_index=torch.tensor(edge_index, dtype=torch.long).t().contiguous(),
                        edge_attr=torch.tensor(edge_attr, dtype=torch.float),
                        edge_type=torch.tensor(edge_type, dtype=torch.long),
                        g=torch.tensor([[1]], dtype=torch.float))
            data_list.append(data)
        return data_list

    @classmethod
    def get_data_list_for_bde_and_map(cls, rdmols, smileses):
        mid_to_num_data = {}
        data_lists_for_bde = []
        for i, rdmol in enumerate(rdmols):
            if rdmol is None:
                mid_to_num_data[i] = -1
                with open(BTPath.CBE_NON_RD_SMILES_FP, 'a', encoding='utf-8')as f:
                    f.write(f"{smileses[i]}\n")
                continue
            try:
                data_list_for_bde = cls.get_data_list_to_bde(rdmol)
                mid_to_num_data[i] = len(data_list_for_bde)
                data_lists_for_bde.extend(data_list_for_bde)
            except Exception as e:
                # raise e
                print(e)
                with open(BTPath.CBE_NON_BDE_SMILES_FP, 'a', encoding='utf-8')as f:
                    f.write(f"{smileses[i]}\n")
                mid_to_num_data[i] = 0

        return data_lists_for_bde, mid_to_num_data

    @classmethod
    def _get_bdes_from_data_list(cls, data_list, num):
        n = 0
        res = None
        while True:
            n += 1
            start_n = (n-1)*num
            end_n = n * num
            if end_n >= len(data_list):
                data = data_list[start_n:]
            else:
                data = data_list[start_n:end_n]
            if len(data) == 0:
                break
            batch_data = Batch.from_data_list(data).to(cls._device)
            bdes = cls._bond_energy_model(batch_data, len(data))
            if res is None:
                res = bdes
            else:
                res = torch.cat((res, bdes), 0)
            if end_n >= len(data_list):
                break
        return res

    @classmethod
    def get_bdes_from_data_list(cls, data_list):
        for num in [6000, 4000, 2000, 1000, 500, 100]:
            try:
                return cls._get_bdes_from_data_list(data_list, num)
            except RuntimeError as e:
                print(e)
                num = int(num / 2)

    @classmethod
    def get_rts_data_list_by_rdmols(cls, rdmols, smileses):
        res = []
        data_lists_for_bde, mid_to_num_data = cls.get_data_list_for_bde_and_map(rdmols, smileses)
        print(f"data num: {len(data_lists_for_bde)}")
        # batch_data = Batch.from_data_list(data_lists_for_bde).to(cls._device)
        # bond_energies = cls._bond_energy_model(batch_data, len(data_lists_for_bde))
        bond_energies = cls.get_bdes_from_data_list(data_lists_for_bde)
        tot = 0
        for i, rdmol in enumerate(rdmols):
            num_data = mid_to_num_data[i]
            if num_data == -1:
                continue
            elif num_data == 0:
                rts_data = cls.get_rts_data_by_rdmol_and_100(rdmol)
                rts_data.key = smileses[i]
                res.append(rts_data)
            else:
                # if num_data == 54:
                #     print(num_data)
                bde = bond_energies[tot:num_data+tot]
                tot += num_data
                rts_data = cls.get_rts_data_by_rdmol_and_bes(rdmol, bde.view(1, num_data)[0].tolist())
                rts_data.key = smileses[i]
                res.append(rts_data)
        return res

    @classmethod
    def get_rdmols(cls, smileses: [str]):
        res = []
        for smiles in smileses:
            try:
                rdmol  = AllChem.MolFromSmiles(smiles)
                res.append(rdmol)
            except:
                res.append(None)
        return res

    @classmethod
    def get_data_list_by_smileses(cls, smileses: [str]):
        rdmols = cls.get_rdmols(smileses)
        data_list = cls.get_rts_data_list_by_rdmols(rdmols, smileses)
        return data_list


if __name__ == "__main__":
    BepDataUtils.get_data_by_smiles('OC(C=O)c1cc[nH]c1')
