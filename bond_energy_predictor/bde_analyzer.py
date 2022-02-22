#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Time    : 2021/10/8 15:45
# @Author  : zhangbc0315@outlook.com
# @File    : bde_analyzer.py
# @Software: PyCharm
import os

import pandas as pd
from rdkit.Chem import AllChem
from rdkit.Chem.rdchem import BondType

BDE_TRAIN_FP = 'C:\\Users\\zhang\\OneDrive\\Projects\\RXN-PREDICTOR\\BondEnergy\\train_bonds.csv'
BDE_TEST_FP = 'C:\\Users\\zhang\\OneDrive\\Projects\\RXN-PREDICTOR\\BondEnergy\\test_bonds.csv'


class BdeAnalyzer:

    @classmethod
    def get_bond_type(cls, smiles, connections, edge_types) -> int:
        mol = AllChem.MolFromSmiles(smiles)
        mol = AllChem.AddHs(mol)
        idx = edge_types.index(0)
        aids = connections[idx]
        bond = mol.GetBondBetweenAtoms(aids[0], aids[1])
        return int(bond.GetBondType())

    @classmethod
    def get_bde_from_fp(cls, fp: str):
        df = pd.read_csv(fp)
        for i, row in df.iterrows():
            elements = [row['Element1'], row['Element2']]
            if 'H' in elements:
                continue
            elements = sorted(elements)
            bond_type_int = cls.get_bond_type(row['Molecule'], eval(row['Connections']), eval(row['EdgeTypes']))
            yield f'{elements[0]}-{elements[1]}-{bond_type_int}', row['Energy']

    @classmethod
    def process(cls):
        pass


if __name__ == "__main__":
    pass
