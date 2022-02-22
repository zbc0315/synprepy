#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Time    : 2021/9/3 13:12
# @Author  : zhangbc0315@outlook.com
# @File    : sane_determiner.py
# @Software: PyCharm
import torch
from torch_geometric.data import Batch
from rdkit.Chem import AllChem

from mcts_lite.mcts_path import MctsPath
from nn_utils.gnn_data_utils import GnnDataUtils
from mcts_lite.policy.spectators_predictor import SpectatorsPredictor
from rxn_determiner.gcn_model import GcnModel
from mcts_lite.mcts_config import MctsConfig
from crkitpy.smilesparse.rdkit_smiles_parser import RdkitSmilesParser
from crkitpy.smatesparse.smates_parser import SmatesParser
from crkitpy.inchiparse.inchi_writer import InchiWriter
from crkitpy.rxnhandler.rxn_template_apply import RxnTemplateApply


class SaneDeterminer:
    cats_len = 325
    sols_len = 321

    device = MctsConfig.ml_device

    model = GcnModel(6, 12, 650).to(device)
    state = torch.load(MctsPath.MODEL_RD_FP, map_location=torch.device(device))
    model.load_state_dict(state['model'])
    model.eval()

    model_without_cs = GcnModel(6, 12, 0).to(device)
    state_without_cs = torch.load(MctsPath.MODEL_RD_WITHOUT_CS_FP, map_location=torch.device(device))
    model_without_cs.load_state_dict(state_without_cs['model'])
    model_without_cs.eval()

    sane_score = 0.45
    # sane_score = 0.5

    # region ===== has other products =====

    @classmethod
    def is_rxn_sane_without_condition(cls, rd_reactants, p_inchi):
        r_data = GnnDataUtils.gnn_data_from_mols(rd_reactants, '')
        r_data.g = torch.tensor([], dtype=torch.float)
        p_data = GnnDataUtils.gnn_data_from_inchi(p_inchi)
        r_batch_data = Batch.from_data_list([r_data]).to(cls.device)
        p_batch_data = Batch.from_data_list([p_data]).to(cls.device)
        pred = cls.model_without_cs(r_batch_data, p_batch_data)
        return pred

    @classmethod
    def _get_rxn_smiles(cls, rd_reactants, rd_product):
        reactants_smiles = '.'.join([AllChem.MolToSmiles(mol) for mol in rd_reactants])
        product_smiles = AllChem.MolToSmiles(rd_product)
        return f'{reactants_smiles}>>{product_smiles}'

    @classmethod
    def _get_crmols_inchis(cls, crmols):
        res = []
        for crmol in crmols:
            try:
                res.append(InchiWriter.mol_to_inchi(crmol))
            except Exception as e:
                print(e)
                return []
        return res

    @classmethod
    def _has_other_products(cls, rd_reactants, rd_product, product_inchi, smates, tid):
        return False
        reactants = [RdkitSmilesParser.rdmol_to_mol(rd_mol) for rd_mol in rd_reactants]
        for products in RxnTemplateApply.reactants_to_products(smates, reactants, strict=True):
            product_inchis = cls._get_crmols_inchis(products)
            if len(product_inchis) > 0 and product_inchi not in product_inchis:
                print(f'has other products: {cls._get_rxn_smiles(rd_reactants, rd_product)},\n{product_inchis}, \n{tid}: {smates}')
                # for p in RxnTemplateApply.reactants_to_products(smates, reactants, strict=True):
                #     print('*'*1000)
                return True
        return False

    # endregion

    @classmethod
    def _get_ids(cls, num):
        num_str = str(num)
        if len(num_str) == 1:
            return 0, num
        elif len(num_str) == 2:
            return int(num_str[0]), int(num_str[1])
        else:
            raise ValueError('预测的催化剂和溶剂过多')

    @classmethod
    def _cats_sols_one_hot_iter(cls, cats, sols):
        """ cat one hot,
                idx = 0 : 无催化剂
                idx = cid+1 : cid催化剂
                idx = len(cats) + 1 : 未知催化剂

        """
        num = (len(cats) + 2) * (len(sols) + 2)

        for n in range(num):
            id1, id2 = cls._get_ids(n)
            if id1 < len(cats):
                cat_oh_idx = cats[id1] + 1
            else:
                cat_oh_idx = cls.cats_len + 1 if id1 == len(cats) else 0
            if id2 < len(sols):
                sol_oh_idx = sols[id2] + 1
            else:
                sol_oh_idx = cls.sols_len + 1 if id2 == len(cats) else 0
            cat_one_hot = [0] * (cls.cats_len + 2)
            sol_one_hot = [0] * (cls.sols_len + 2)
            cat_one_hot[cat_oh_idx] = 1
            sol_one_hot[sol_oh_idx] = 1
            yield torch.tensor(cat_one_hot, dtype=torch.float), \
                  torch.tensor(sol_one_hot, dtype=torch.float), \
                  cat_oh_idx, sol_oh_idx

    @classmethod
    def get_r_gnn_data(cls, rd_reactants, cats, sols):
        r_data = GnnDataUtils.gnn_data_from_mols(rd_reactants, '')
        r_data_list = []
        cat_sol_oh_idx_list = []
        for cat_one_hot, sol_one_hot, cat_id, sol_id in cls._cats_sols_one_hot_iter(cats, sols):
            global_vector = torch.cat([cat_one_hot, sol_one_hot])
            global_vector = torch.reshape(global_vector, [1, len(cat_one_hot) + len(sol_one_hot)])
            res = r_data.clone()
            res.g = global_vector
            r_data_list.append(res)
            cat_sol_oh_idx_list.append((cat_id, sol_id))
        return r_data_list, cat_sol_oh_idx_list

    @classmethod
    def _get_p_batch_data(cls, p_data, num):
        p_data_list = [p_data.clone() for i in range(num)]
        return Batch.from_data_list(p_data_list)

    @classmethod
    def _get_sane_reaction_cat_sol(cls, pred, cat_sol_list):
        res = []
        sane_idx = torch.where(pred >= cls.sane_score)[0].tolist()
        for sane_id in sane_idx:
            res.append(cat_sol_list[sane_id])
        return res

    @classmethod
    def is_sane_consider_spec(cls, rd_reactants_list, rd_product, tids, smateses):
        product_inchi = AllChem.MolToInchi(rd_product)
        cats_list, sols_list = SpectatorsPredictor.batch_predict_cats_sols(rd_reactants_list, rd_product, tids)
        p_data = GnnDataUtils.gnn_data_from_mol(rd_product)
        for n, rd_reactants in enumerate(rd_reactants_list):
            tid = tids[n]
            smates = smateses[n]
            cats = cats_list[n]
            sols = sols_list[n]
            r_data_list, cat_sol_oh_idx_list = cls.get_r_gnn_data(rd_reactants, cats, sols)
            r_data_batch = Batch.from_data_list(r_data_list).to(cls.device)
            p_data_batch = cls._get_p_batch_data(p_data, len(r_data_list)).to(cls.device)
            pred = cls.model(r_data_batch, p_data_batch)

            sane_cat_sol_oh_idx_list = cls._get_sane_reaction_cat_sol(pred, cat_sol_oh_idx_list)
            if len(sane_cat_sol_oh_idx_list) == 0 or cls._has_other_products(rd_reactants, rd_product, product_inchi, smates, tid):
                yield False, [], tid, [], []
            else:
                yield True, sane_cat_sol_oh_idx_list, tid, cats, sols


def test():
    r11 = 'COc1cncc2[nH]ccc12'
    r12 = 'O=C1CCC(=O)N1Cl'

    r21 = 'CO'
    r22 = 'CN'

    r11 = 'O=C(C=C1)N(C2=CC=CC=C2)C1=O'
    r12 = 'BrCC(C(OCC)=O)=CC'

    p = 'O=C1N(C2=CC=CC=C2)C(C3CC=C(C(OCC)=O)CC31)=O'
    # p = 'COc1cncc2[nH]cc(Cl)c12'

    rd11 = AllChem.MolFromSmiles(r11)
    rd12 = AllChem.MolFromSmiles(r12)

    rd21 = AllChem.MolFromSmiles(r21)
    rd22 = AllChem.MolFromSmiles(r22)

    prd = AllChem.MolFromSmiles(p)

    for is_sane, cat_sol_list, tid in SaneDeterminer.is_sane_consider_spec([[rd11, rd12]], prd, [100]):
        print(cat_sol_list)
        print(f'{is_sane}: {len(cat_sol_list)}')


def test1():
    rs_smileses = ['CCOC(C(c(cn1COP(O)(O)=O)c2c1c(n3nc(C)nc3)ncc2OC)=O)=O', 'O=C(N1CCNCC1)c2ccccc2']
    p_inchi = 'InChI=1S/C13H18N2O/c1-2-14-8-10-15(11-9-14)13(16)12-6-4-3-5-7-12/h3-7H,2,8-11H2,1H3'
    rs_mols = [AllChem.MolFromSmiles(smiles) for smiles in rs_smileses]
    print(SaneDeterminer.is_rxn_sane_without_condition(rs_mols, p_inchi))


if __name__ == "__main__":
    test1()
