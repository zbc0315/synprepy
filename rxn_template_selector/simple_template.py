# -*- coding: utf-8 -*-
# @Time     : 2021/4/24 21:20
# @Author   : zbc
# @Email    : zhangbc0315@outlook.com
# @FileName : simple_template.py
# @Software : PyCharm

from rdkit.Chem import AllChem
from rdkit.Chem.AllChem import ChemicalReaction

from rxn_template_selector.rxn_center_detector import RxnCenterDetector
from rxn_template_selector.mol_handler import MolHandler


class SimpleTemplate:

    # =============
    # -- extract --

    @classmethod
    def load_rxn(cls, rxn_smiles: str):
        reactant = AllChem.MolFromSmiles(rxn_smiles.split('>')[0])
        product = AllChem.MolFromSmiles(rxn_smiles.split('>')[-1])
        rxn = AllChem.ChemicalReaction()
        rxn.AddReactantTemplate(reactant)
        rxn.AddProductTemplate(product)
        return rxn

    @classmethod
    def extract_by_smiles(cls, rxn_smiles) -> str:
        # return cls.extract_by_rxn(cls.load_rxn(rxn_smiles))
        return cls.extract_by_rxn(AllChem.ReactionFromSmarts(rxn_smiles))

    @classmethod
    def extract_by_rxn(cls, rxn) -> str:
        reactant_center, product_center = RxnCenterDetector.detect_rxn_center(rxn)
        reactants = list(rxn.GetReactants())
        products = list(rxn.GetProducts())
        reactants_smarts = [MolHandler.smarts_of_frag_by_aids(reactants[r], reactant_center[r])
                            for r in reactant_center.keys()]
        products_smarts = [MolHandler.smarts_of_frag_by_aids(products[p], product_center[p])
                           for p in product_center.keys()]
        return f"{'.'.join(reactants_smarts)}>>{'.'.join(products_smarts)}"

    # ==================================
    # -- apply: reactants to products --

    @staticmethod
    def _load_template(temp_smarts: str):
        r = AllChem.MolFromSmarts(temp_smarts.split('>>')[0])
        p = AllChem.MolFromSmarts(temp_smarts.split('>>')[-1])
        rxn = ChemicalReaction()
        rxn.AddReactantTemplate(r)
        rxn.AddProductTemplate(p)
        temp_smarts = temp_smarts.split('>')
        temp_smarts.reverse()
        temp_smarts = '>'.join(temp_smarts)
        return AllChem.ReactionFromSmarts(temp_smarts)

    @classmethod
    def reactants_to_products_by_mols_to_mols(cls, rxn_temp, reactants):
        for products in rxn_temp.RunReactants(reactants):
            yield products

    @classmethod
    def reactants_to_products_by_smarts_to_mols(cls, temp_smarts: str, reactants_smarts: [str]):
        rxn_temp = cls._load_template(temp_smarts)
        reactants = [AllChem.MolFromSmarts(smarts) for smarts in reactants_smarts]
        for products in cls.reactants_to_products_by_mols_to_mols(rxn_temp, reactants):
            yield products

    @classmethod
    def reactants_to_products_by_smarts_to_smarts(cls, temp_smarts: str, reactants_smarts: [str]):
        for products in cls.reactants_to_products_by_smarts_to_mols(temp_smarts, reactants_smarts):
            yield [AllChem.MolToSmiles(product) for product in products]


if __name__ == '__main__':
    smi = "[CH:1]1=[CH:6][CH:5]=[CH:4][CH:3]=[CH:2]1.[CH3:8][Cl:7]>>[CH3:8][C:1]1=[CH:6][CH:5]=[CH:4][CH:3]=[CH:2]1"
    # print()
    sma = SimpleTemplate.extract_by_smiles(smi)
    print(sma)
    # rs_smi = ["CCCCCCl"]
    rs_smi = ["CCC1=CC(=CC(C)=C1)C1CC1"]
    for ps in SimpleTemplate.reactants_to_products_by_smarts_to_smarts(sma, rs_smi):
        for p in ps:
            print(p)
