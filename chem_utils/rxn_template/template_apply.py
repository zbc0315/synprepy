#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Time    : 2022/1/5 11:01
# @Author  : zhangbc0315@outlook.com
# @File    : template_apply.py
# @Software: PyCharm

from rdkit.Chem.rdchem import Mol
from rdkit.Chem import AllChem


class TemplateApply:

    @classmethod
    def synthesis_result_to_smiles(cls, syn_results: [[Mol]]) -> [[str]]:
        res = []
        for mols in syn_results:
            smileses = []
            for mol in mols:
                AllChem.SanitizeMol(mol)
                smileses.append(AllChem.MolToSmiles(mol))
            if smileses not in res:
                res.append(smileses)
        return res

    @classmethod
    def synthesis(cls, rxn_temp: AllChem.ChemicalReaction, reactant: Mol):
        res = rxn_temp.RunReactants((reactant,))
        for ms in res:
            for m in ms:
                AllChem.RemoveStereochemistry(m)
        return res

    @classmethod
    def synthesis_by_smarts(cls, rxn_temp_smarts: str, reactant_smiles: str):
        return cls.synthesis(AllChem.ReactionFromSmarts(rxn_temp_smarts), AllChem.MolFromSmiles(reactant_smiles))

    @classmethod
    def retro_synthesis_by_smarts(cls, rxn_temp_smarts: str, product_smiles: str):
        retro_rxn_temp_smarts = f"{rxn_temp_smarts.split('>>')[-1]}>>{rxn_temp_smarts.split('>>')[0]}"
        return cls.synthesis_by_smarts(retro_rxn_temp_smarts, product_smiles)

    @classmethod
    def retro_synthesis_by_smarts_to_smiles(cls, rxn_temp_smarts: str, product_smiles: str):
        reactants = cls.retro_synthesis_by_smarts(rxn_temp_smarts, product_smiles)
        return [[AllChem.MolToSmiles(m) for m in ms] for ms in reactants]


if __name__ == "__main__":
    p_smiles = "FC1=CC=C(C=C1)CN2C(C)=C(C3=CC(C(NC4CNCC4)=O)=CC=C23)C"
    # p_smiles = "C1CCNC1"
    r_smiles = "CC(C)(OC(N1CC(CC1)NC(C2=CC=C3N(C(C)=C(C3=C2)C)CC4=CC=C(C=C4)F)=O)=O)C.Cl"
    temp_smarts = "(CC(C)(C)OC(=O)[N:1].Cl)>>[NH:1]"
    # temp_smarts = "(CC(C)(C)OC(=O)[N:2]([CH2:1])[CH2:3].Cl)>>[CH2:1][NH:2][CH2:3]"
    # temp_smarts = "C[O:1]>>[OH:1]"
    print(TemplateApply.synthesis_result_to_smiles(TemplateApply.synthesis_by_smarts(temp_smarts, r_smiles)))
    print(TemplateApply.synthesis_result_to_smiles(TemplateApply.retro_synthesis_by_smarts(temp_smarts, p_smiles)))
