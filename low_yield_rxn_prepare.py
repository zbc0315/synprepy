#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Time    : 2022/2/22 11:18
# @Author  : zhangbc0315@outlook.com
# @File    : low_yield_rxn_prepare.py
# @Software: PyCharm

from tqdm import tqdm
import pandas as pd
from rdkit.Chem import AllChem
from rdkit.Chem.rdchem import KekulizeException
from rdkit import RDLogger

RDLogger.DisableLog('rdApp.*')

from chem_utils.chem_handler import MolHandler
from chem_utils.rxn_template.template_apply import TemplateApply
from config import RDConfig, Config
from data_utils import RxnData, RxnDataType, RxnTemplateData
from utils import TSVUtils


class LowYieldRxnPrepare:

    def __init__(self, rd_config: RDConfig):
        self._config = rd_config
        self._rxn_data = RxnData(rd_config.rxn_data_fp, None, RxnDataType.TRAIN_TEST, rd_config.evaluate_year)
        self._rxn_data.init_yield(95, 100)
        self._template_data = RxnTemplateData(rd_config.rct_fp)
        self._template_data.init_filter(100)

    @classmethod
    def _get_main_product(cls, ps_inchis: [str]):
        res = ps_inchis[0]
        if len(ps_inchis) == 1:
            return res
        else:
            for p_inchi in ps_inchis:
                if len(p_inchi) > len(res):
                    res = p_inchi
            return res

    @classmethod
    def _get_rxn_smiles(cls, reactant_smiles: str, product_inchi: str):
        r_mol = AllChem.MolFromSmiles(reactant_smiles)
        AllChem.RemoveStereochemistry(r_mol)
        r_smiles = MolHandler.get_smiles_without_map_num(r_mol)
        p_smiles = AllChem.MolToSmiles(AllChem.MolFromInchi(product_inchi))
        return r_smiles + ">>" + p_smiles

    def _get_diff_products(self, reactant_smiles, product_smiles):
        p_mol = AllChem.MolFromSmiles(product_smiles)
        AllChem.RemoveStereochemistry(p_mol)
        product_inchi = AllChem.MolToInchi(p_mol)
        for rt in self._template_data.get_all_template():
            res = TemplateApply.synthesis_by_smarts(rt.rt_smiles, reactant_smiles)
            for ps in res:
                try:
                    ps_inchis = [AllChem.MolToInchi(p) for p in ps]
                except KekulizeException as e:
                    continue
                if product_inchi in ps_inchis:
                    continue
                main_p_inchi = self._get_main_product(ps_inchis)
                yield main_p_inchi, rt.rt_code

    def process(self):
        df = pd.DataFrame(columns=["ly_rxn_code", "rxn_code", "rt_code", "rxn_smiles", "catalysts_codes", "solvents_codes"])
        i = 0
        num = 0
        with tqdm(total=len(self._rxn_data))as pbar:
            for rxn in self._rxn_data.get_all_rxn(False):
                pbar.set_postfix_str(f"get {num} rxns")
                pbar.update(1)
                reactant_smiles, _, product_smiles = rxn.rxn_smiles.split('>')
                diff_ps = set()
                for diff_p_inchi, rt_code in self._get_diff_products(reactant_smiles, product_smiles):
                    if diff_p_inchi in diff_ps:
                        continue
                    diff_ps.add(diff_p_inchi)
                    try:
                        rxn_smiles = self._get_rxn_smiles(reactant_smiles, diff_p_inchi)
                    except:
                        continue
                    num += 1
                    df = df.append({"ly_rxn_code": f"LR{i}",
                                    "rxn_code": rxn.rxn_code,
                                    "rt_code": rt_code,
                                    "rxn_smiles": rxn_smiles,
                                    "catalysts_codes": rxn.catalysts_codes,
                                    "solvents_codes": rxn.solvents_codes}, ignore_index=True)
                    i += 1
                    if len(df) >= 1000:
                        TSVUtils.df_to_tsv(df, self._config.low_yield_rxns_tsv_fp, 'a')
                        df = pd.DataFrame(columns=["ly_rxn_code", "rxn_code", "rt_code", "rxn_smiles", "catalysts_codes",
                                                   "solvents_codes"])
            if len(df) > 0:
                TSVUtils.df_to_tsv(df, self._config.low_yield_rxns_tsv_fp, 'a')


if __name__ == "__main__":
    LowYieldRxnPrepare(Config("config.json").rd_config).process()
