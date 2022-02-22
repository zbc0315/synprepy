#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Time    : 2022/2/22 10:54
# @Author  : zhangbc0315@outlook.com
# @File    : by_products_predictor.py
# @Software: PyCharm

import torch
from torch_geometric.data import Batch
from rdkit.Chem import AllChem
from rdkit import RDLogger

RDLogger.DisableLog('rdApp.*')

from config import Config, RCTSConfig
from utils.gnn_data_utils import GNNDataUtils
from data_utils import SpectatorData
from utils import ModelUtils
from rxn_template_selector.rct_predictor import RCTPredictor
from chem_utils.rxn_template.template_apply import TemplateApply


class ByProductsPredictor:

    def __init__(self, rcts_config: RCTSConfig):
        self._cat_data = SpectatorData(rcts_config.cat_code_and_ohi_fp)
        self._sol_data = SpectatorData(rcts_config.sol_code_and_ohi_fp)
        self._predictor = RCTPredictor(rcts_config)

    def predict(self, reactant_smiles: str, product_smiles: str, cat_code: str, sol_code: str):
        p_mol = AllChem.MolFromSmiles(product_smiles)
        AllChem.RemoveStereochemistry(p_mol)
        p_inchi = AllChem.MolToInchi(p_mol)
        data = GNNDataUtils.get_gnn_data_from_smiles(reactant_smiles)
        data.g = torch.tensor([self._cat_data.get_one_hot_vec_by_mol_code(cat_code) +
                               self._sol_data.get_one_hot_vec_by_mol_code(sol_code)], dtype=torch.float)
        batch = Batch.from_data_list([data])
        rt_smiles_list = self._predictor.predict_to_rt_smiles(batch, 5)
        predicted_by_ps = set()
        for rt_smiles in rt_smiles_list:
            for products in TemplateApply.synthesis_by_smarts(rt_smiles, reactant_smiles):
                for product in products:
                    AllChem.RemoveStereochemistry(product)
                by_ps_inchi = [AllChem.MolToInchi(product) for product in products]
                if p_inchi in by_ps_inchi:
                    continue
                if any([bp_inchi in predicted_by_ps for bp_inchi in by_ps_inchi]):
                    continue
                else:
                    predicted_by_ps = predicted_by_ps.union(set(by_ps_inchi))
                    print("\n".join(by_ps_inchi))
                    print("*"*100)


if __name__ == "__main__":
    r = "O=C(c1ccccc1)NC(N2OCC(CN(C3)C(OCc4ccccc4)=O)C23c5cccc(Br)c5)=S"
    p = "O=C(NC1=NC2(C(CN(C2)C(OCc3ccccc3)=O)CS1)c4cccc(Br)c4)c5ccccc5"
    r = "O=C1CCC(N1Br)=O.O=C2c3c(CCC2)cc(C#N)cc3"
    p = "O=C1c2c(CCC1Br)cc(C#N)cc2"
    r = "OB(c1ccc(c2ccccc2)cc1)O.CC(C)(OC(c3c(c4ccc(F)cc4)c(CNCC5(C)C)c5nc3C(C)C)C(OC)=O)C"
    p = "CC(C)(OC(c1c(c2ccc(F)cc2)c(CN(CC3(C)C)c4ccc(c5ccccc5)cc4)c3nc1C(C)C)C(OC)=O)C"
    r = "CC(OC(C)=O)=O.CC(C)(OC(c1c(c2ccc(F)cc2)c(CNCC3(C)C)c3nc1C(C)C)C(O)=O)C"
    p = "CC(C)(OC(c1c(c2ccc(F)cc2)c(CN(CC3(C)C)C(C)=O)c3nc1C(C)C)C(O)=O)C"
    c = None
    s = "M821"
    r = "BrC(C)(C)C(OCC)=O.Sc1cccnc1c2ccc(C#N)cc2"
    p = "N#Cc1ccc(c2ncccc2SC(C)(C)C(OCC)=O)cc1"
    c = None
    s = "M47"
    r = "Brc1cncnc1c(ccc2C#N)c3c2cccc3.SCC(OC)=O"
    p = "O=C(OC)CSc1cncnc1c(ccc2C#N)c3c2cccc3"
    c = None
    s = "M47"
    r = "Clc(c(C(OC)=O)c(=O)[nH]1)c2c1cccc2.COc3ccc(CN)cc3"
    p = "COc1ccc(CNc(c(C(OC)=O)c(=O)[nH]2)c3c2cccc3)cc1"
    c = None
    s = "M47"
    ByProductsPredictor(Config("../config.json").rcts_config).predict(r, p, c, s)
