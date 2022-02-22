#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Time    : 2021/9/3 10:13
# @Author  : zhangbc0315@outlook.com
# @File    : check_reactants.py
# @Software: PyCharm

from rdkit.Chem import AllChem

from crkitpy.inchiparse.inchi_writer import InchiWriter
from mcts_lite.mcts_error import *
from log_utils.logger import Logging


class CheckReactants:

    @classmethod
    def _mols_to_inchis(cls, mols):
        inchis = [InchiWriter.mol_to_inchi(mol) for mol in mols]
        return inchis

    @classmethod
    def _inchis_to_rdmols(cls, inchis):
        rdmols = [AllChem.MolFromInchi(inchi) for inchi in inchis]
        if any([rdmol is None for rdmol in rdmols]):
            raise UnexpectedError(f'inchi 无法解析为rdmol')
        return rdmols

    @classmethod
    def _check_reactants_not_contain_product(cls, rs_inchis, p_inchi):
        if p_inchi in rs_inchis:
            raise UnexpectedError(f'逆反应中，反应物包含了产物')

    @classmethod
    def _check_reactants_not_contain_target(cls, rs_inchis, target_inchi):
        if target_inchi in rs_inchis:
            raise ExpectedError(f'逆反应中，反应物包含了目标化合物')

    @classmethod
    def _check_reactants_not_yielded(cls, rs_inchis, yielded_reactants):
        if set(rs_inchis) in yielded_reactants:
            raise ExpectedError(f'逆反应中，反应物集合已经被采纳至少一次了')

    @classmethod
    def _check_reactants_in_intermediate_products(cls, rs_inchis, past_products_inchis):
        if any([r_inchi in past_products_inchis for r_inchi in rs_inchis]):
            raise ExpectedError(f'你反应中，反应物包含了前述反应的产物')

    @classmethod
    def _check_special_reactants(cls, reactants_inchis):
        wrong_reactants = [
                            'InChI=1S/CH4/h1H4',  # CH4
                            'InChI=1S/N2/c1-2',  # N2
                            'InChI=1S/O2/c1-2',  # O2
                            'InChI=1S/F2/c1-2',  # F2
                            'InChI=1S/C2Cl2O2/c3-1(5)2(4)6',
                            'InChI=1S/BrH/h1H',
                           ]
        if any([inchi in wrong_reactants for inchi in reactants_inchis]):
            raise ExpectedError(f'反应物中，包含不合适的反应物，例如常温下的气态物质，不稳定的物质')

    @classmethod
    def check(cls, reactants, product_inchi, target_inchi, yielded_reactants, intermediate_products=None):
        try:
            reactants_inchis = cls._mols_to_inchis(reactants)
            cls._check_special_reactants(reactants_inchis)
            rd_reactants = cls._inchis_to_rdmols(reactants_inchis)
            cls._check_reactants_not_contain_target(reactants_inchis, target_inchi)
            cls._check_reactants_not_contain_product(reactants_inchis, product_inchi)
            cls._check_reactants_not_yielded(reactants_inchis, yielded_reactants)
            cls._check_reactants_in_intermediate_products(reactants_inchis, intermediate_products)
            return reactants_inchis, rd_reactants, True
        except ExpectedError as e:
            Logging.log.debug(e)
            return None, None, False
        except Exception as e:
            Logging.log.info(e)
            return None, None, False


if __name__ == "__main__":
    pass
