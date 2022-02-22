# -*- coding: utf-8 -*-
# @Time     : 2021/4/27 10:20
# @Author   : zbc
# @Email    : zhangbc0315@outlook.com
# @FileName : gnn_data_utils.py
# @Software : PyCharm

import os

import pandas as pd
import torch
from torch_geometric import data as gdata
from rdkit.Chem import AllChem

from rdkit import RDLogger

RDLogger.DisableLog('rdApp.*')


class GnnDataUtils:

    _periodic_df = None

    @staticmethod
    def _int_from_bool(b):
        """ True = 1, False = 0

        :param b:
        :return:
        """
        return 1 if b else 0

    # @staticmethod
    # def _ont_hot_bond_type(bond_type):
    #     """ 将键的类型转化为 one hot 向量
    #
    #     :param bond_type:
    #     :return:
    #     """
    #     res = [0]*13
    #     res[int(bond_type)] = 1
    #     return res

    @staticmethod
    def _one_hot_bond_type(bond_type):
        bond_level = int(bond_type)
        if bond_level == 12:
            bond_level = 2
        elif bond_level >= 2:
            bond_level += 1
        res = [0]*12
        res[bond_level] = 1
        return res

    # ===============================

    @classmethod
    def bond_to_h(cls, bond):
        return bond.GetBeginAtom().GetAtomicNum() == 1 or bond.GetEndAtom().GetAtomicNum() == 1

    @classmethod
    def change_aid(cls, aid, old_aid_to_new_aid):
        if old_aid_to_new_aid is None:
            return aid
        else:
            return old_aid_to_new_aid[aid]

    @classmethod
    def connect_info_from_mol(cls, mol, no_h=False, old_aid_to_new_aid=None):
        """ 获得键连信息

        :param no_h:
        :param old_aid_to_new_aid:
        :param mol:
        :return:
        """
        connect_info = []
        for bond in mol.GetBonds():
            if no_h and cls.bond_to_h(bond):
                continue
            baid = bond.GetBeginAtomIdx()
            eaid = bond.GetEndAtomIdx()
            connect_info.append([cls.change_aid(baid, old_aid_to_new_aid),
                                 cls.change_aid(eaid, old_aid_to_new_aid)])
            connect_info.append([cls.change_aid(eaid, old_aid_to_new_aid),
                                 cls.change_aid(baid, old_aid_to_new_aid)])
        return connect_info

    # =====================================

    @classmethod
    def _edge_features_from_bond(cls, bond, predict_type) -> [int]:
        """ 获得键的特征，是键的类型的 one hot 向量
        TODO 此处可用 cls._ont_hot_bond_type: 预测催化剂、溶剂、考虑键能
               也可用 cls._one_hot_bond_type: highway_gcn

        :param bond:
        :return:
        """
        return cls._one_hot_bond_type(bond.GetBondType())
        # if predict_type == 'cat_sol':
        #     return cls._ont_hot_bond_type(bond.GetBondType())
        # else:
        #     return cls._one_hot_bond_type(bond.GetBondType())

    @classmethod
    def edges_features_from_mol(cls, mol, predict_type, no_h=False):
        """ 获得所有键的特征

        :param predict_type:
        :param no_h:
        :param mol:
        :return:
        """
        if not no_h:
            features: [[int]] = []
            for bond in mol.GetBonds():
                features.append(cls._edge_features_from_bond(bond, predict_type))
                features.append(cls._edge_features_from_bond(bond, predict_type))
            return features
        else:
            features: [[int]] = []
            old_bid_to_new_bid = {}
            now_bid = 0
            for i, bond in enumerate(mol.GetBonds()):
                if cls.bond_to_h(bond):
                    continue
                features.append(cls._edge_features_from_bond(bond, predict_type))
                features.append(cls._edge_features_from_bond(bond, predict_type))
                old_bid_to_new_bid[i] = now_bid
                now_bid += 1
            return features, old_bid_to_new_bid

    # =========================================

    @classmethod
    def _load_periodic_data(cls):
        fp = '../../static/periodic.tsv'
        if not os.path.exists(fp):
            fp = '../static/periodic.tsv'
        df = pd.read_csv(fp, sep='\t', index_col=0)
        return df

    @classmethod
    def _hybrid_to_num(cls, hybrid):
        id_hybrid = int(hybrid)
        return id_hybrid/10

    @classmethod
    def _bool_to_num(cls, b: bool) -> int:
        return 1 if b else 0

    @classmethod
    def _node_features_from_atom(cls, atom, predict_type, transform_ex_h_to_im_h: bool = False) -> [int]:
        """ 获得原子的特征
        TODO 此处可加 atom.GetFormalCharge: 预测催化剂、溶剂、考虑键能
             也可不加 atom.GetFormalCharge: highway_gcn

        :param atom:
        :return:
        """
        if cls._periodic_df is None:
            cls._periodic_df = cls._load_periodic_data()
        if transform_ex_h_to_im_h:
            hs_num = 0
            for nbr_atom in atom.GetNeighbors():
                if nbr_atom.GetAtomicNum() == 1:
                    hs_num += 1
        else:
            hs_num = atom.GetTotalNumHs()
        return [
            atom.GetAtomicNum(),
            hs_num,
            atom.GetFormalCharge(),
            atom.GetNumRadicalElectrons(),
            atom.GetTotalValence(),
            cls._hybrid_to_num(atom.GetHybridization()),
            cls._bool_to_num(atom.IsInRing()),
            cls._bool_to_num(atom.GetIsAromatic()),
            cls._periodic_df.loc[atom.GetAtomicNum(), 'row'],
            cls._periodic_df.loc[atom.GetAtomicNum(), 'column'],
        ]
        # if predict_type == 'cat_sol':
        #     return [atom.GetAtomicNum(),
        #             atom.GetFormalCharge(),
        #             hs_num,
        #             atom.GetTotalValence(),
        #             atom.GetNumRadicalElectrons(),
        #             cls._int_from_bool(atom.IsInRing()),
        #             cls._int_from_bool(atom.GetIsAromatic())]
        # else:
        #     return [atom.GetAtomicNum(),
        #             # atom.GetFormalCharge(),
        #             hs_num,
        #             atom.GetTotalValence(),
        #             atom.GetNumRadicalElectrons(),
        #             cls._int_from_bool(atom.IsInRing()),
        #             cls._int_from_bool(atom.GetIsAromatic())]

    @classmethod
    def nodes_features_from_mol(cls, mol, predict_type: str, no_h: bool = False,
                                transform_ex_h_to_im_h: bool = False) -> [[int]]:
        """ 获得所有原子的特征

        :param predict_type:
        :param no_h:
        :param mol:
        :return:
        """
        if not no_h:
            features: [(int, [int])] = []
            for atom in mol.GetAtoms():
                features.append((atom.GetIdx(), cls._node_features_from_atom(atom, predict_type, transform_ex_h_to_im_h)))
            sorted_features = sorted(features, key=lambda x: x[0])
            return [sorted_feature[1] for sorted_feature in sorted_features]
        else:
            features = []
            old_aid_to_new_aid = {}
            now_aid = 0
            for atom in mol.GetAtoms():
                if atom.GetAtomicNum() == 1:
                    continue
                old_aid_to_new_aid[atom.GetIdx()] = now_aid
                now_aid += 1
                features.append(cls._node_features_from_atom(atom, predict_type, transform_ex_h_to_im_h))
            return features, old_aid_to_new_aid

    @classmethod
    def gnn_data_from_mol(cls, mol, predict_type: str = ''):
        nodes = torch.tensor(cls.nodes_features_from_mol(mol, predict_type), dtype=torch.float)
        edges = torch.tensor(cls.edges_features_from_mol(mol, predict_type), dtype=torch.float)
        connect = torch.tensor(cls.connect_info_from_mol(mol), dtype=torch.long)
        data = gdata.Data(x=nodes, edge_index=connect.t().contiguous(), edge_attr=edges)
        return data

    @classmethod
    def gnn_data_from_mols(cls, mols, predict_type):
        mol = mols[0]
        for n, another_mol in enumerate(mols):
            if n == 0:
                continue
            mol = AllChem.CombineMols(mol, another_mol)
        return cls.gnn_data_from_mol(mol, predict_type)

    @classmethod
    def gnn_data_from_smiles(cls, smiles: str):
        mol = AllChem.MolFromSmiles(smiles)
        if mol is None:
            return None
        else:
            return cls.gnn_data_from_mol(mol)

    @classmethod
    def _is_radical_mol(cls, mol):
        for atom in mol.GetAtoms():
            if atom.GetNumRadicalElectrons() > 0:
                return True
        return False

    @classmethod
    def gnn_data_from_inchi(cls, inchi: str, no_radical: bool = False):
        try:
            mol = AllChem.MolFromInchi(inchi)
        except:
            return None
        if mol is None:
            return None
        elif no_radical and cls._is_radical_mol(mol):
            return None
        else:
            return cls.gnn_data_from_mol(mol)

    @classmethod
    def gnn_data_from_smileses(cls, smileses: [str]):
        smiles = '.'.join(smileses)
        return cls.gnn_data_from_smiles(smiles)

    @classmethod
    def gnn_data_from_inchis(cls, inchis: str):
        try:
            smileses = [AllChem.MolToSmiles(AllChem.MolFromInchi(inchi)) for inchi in inchis]
        except:
            return None
        return cls.gnn_data_from_smileses(smileses)


if __name__ == '__main__':
    pass
