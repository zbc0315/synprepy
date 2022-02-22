#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Time    : 2022/1/8 21:26
# @Author  : zhangbc0315@outlook.com
# @File    : gnn_data_utils.py
# @Software: PyCharm

import torch
from torch_geometric.data import Data as GData
from rdkit.Chem.rdchem import Atom, Bond, Mol
from rdkit.Chem import AllChem

from utils.periodic import Periodic


class GNNDataUtils:

    _periodic = Periodic()
    _kernels = [86, 54, 36, 18, 10, 2]

    # region ===== node =====

    @classmethod
    def _bool_to_num(cls, b: bool) -> int:
        return 1 if b else 0

    @classmethod
    def _hybrid_to_num(cls, hybrid) -> int:
        id_hybrid = int(hybrid)
        return id_hybrid

    @classmethod
    def _get_node_features_from_rdatom(cls, atom: Atom) -> [int]:
        row, col = cls._periodic.get_row_and_col(atom.GetAtomicNum())
        return [atom.GetAtomicNum(),
                atom.GetTotalNumHs(),
                atom.GetFormalCharge(),
                atom.GetNumRadicalElectrons(),
                atom.GetTotalValence(),
                cls._hybrid_to_num(atom.GetHybridization()),
                cls._bool_to_num(atom.IsInRing()),
                cls._bool_to_num(atom.GetIsAromatic()),
                row,
                col]

    @classmethod
    def _get_nodes_features_from_rdmol(cls, mol: Mol) -> [[int]]:
        return [cls._get_node_features_from_rdatom(mol.GetAtomWithIdx(i)) for i in range(mol.GetNumAtoms())]

    # endregion

    # region ===== edge & connection =====

    @classmethod
    def _get_edge_features_from_rdbond(cls, bond: Bond) -> [int]:
        bond_level = int(bond.GetBondType())
        if bond_level == 12:
            bond_level = 2
        elif bond_level >= 2:
            bond_level += 1
        res = [0] * 5
        res[bond_level] = 1
        return res

    @classmethod
    def _get_connection_from_rdbond(cls, bond: Bond, reverse: bool = False) -> [int]:
        if reverse:
            return [bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()]
        else:
            return [bond.GetEndAtomIdx(), bond.GetBeginAtomIdx()]

    @classmethod
    def _get_edges_features_and_connections_from_rdmol(cls, mol: Mol) -> ([[int]], [[int]]):
        edges_features = []
        connections = []
        for bond in mol.GetBonds():
            edges_features.append(cls._get_edge_features_from_rdbond(bond))
            edges_features.append(cls._get_edge_features_from_rdbond(bond))
            connections.append(cls._get_connection_from_rdbond(bond, False))
            connections.append(cls._get_connection_from_rdbond(bond, True))
        return edges_features, connections

    # endregion

    @classmethod
    def get_gnn_data_from_rdmol(cls, mol: Mol) -> GData:
        edges_features, connections = cls._get_edges_features_and_connections_from_rdmol(mol)
        x = torch.tensor(cls._get_nodes_features_from_rdmol(mol), dtype=torch.float)
        edge_index = torch.tensor(connections, dtype=torch.long).t().contiguous()
        edge_attr = torch.tensor(edges_features, dtype=torch.float)
        return GData(x=x, edge_index=edge_index, edge_attr=edge_attr)

    @classmethod
    def get_gnn_data_from_smiles(cls, smiles: str) -> GData:
        mol = AllChem.MolFromSmiles(smiles)
        if mol is None:
            print(smiles)
        return cls.get_gnn_data_from_rdmol(mol)

    @classmethod
    def get_gnn_data_from_rxn_smiles(cls, rxn_smiles: str) -> GData:
        r_smiles, _, p_smiles = rxn_smiles.split('>')
        r_mol = AllChem.MolFromSmiles(r_smiles)
        p_mol = AllChem.MolFromSmiles(p_smiles)
        data = GData()
        data.x = torch.tensor(cls._get_nodes_features_from_rdmol(r_mol), dtype=torch.float)
        data.px = torch.tensor(cls._get_nodes_features_from_rdmol(p_mol), dtype=torch.float)
        r_edge_features, r_conn = cls._get_edges_features_and_connections_from_rdmol(r_mol)
        p_edge_features, p_conn = cls._get_edges_features_and_connections_from_rdmol(p_mol)
        data.edge_index = torch.tensor(r_conn, dtype=torch.long).t().contiguous()
        data.p_edge_index = torch.tensor(p_conn, dtype=torch.long).t().contiguous()
        data.edge_attr = torch.tensor(r_edge_features, dtype=torch.float)
        data.p_edge_attr = torch.tensor(p_edge_features, dtype=torch.float)
        return data


if __name__ == "__main__":
    GNNDataUtils.get_gnn_data_from_smiles("CN=C")
