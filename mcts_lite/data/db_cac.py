#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Time    : 2021/8/31 12:43
# @Author  : zhangbc0315@outlook.com
# @File    : db_cac.py
# @Software: PyCharm

import pandas as pd
from rdkit.Chem import AllChem
from rdkit import RDLogger

from mcts_lite.mcts_path import MctsPath

RDLogger.DisableLog('rdApp.*')


class DbCac:

    inchi_fps = {'bbm': MctsPath.CAC_INCHI_NON_STEREO_BBM_FP,
                 'sm': MctsPath.CAC_INCHI_NON_STEREO_SM_FP,
                 'rm': MctsPath.CAC_INCHI_NON_STEREO_RM_FP}

    smiles_fps = {'bbm': MctsPath.CAC_SMILES_NON_STEREO_BBM_FP,
                  'sm': MctsPath.CAC_SMILES_NON_STEREO_SM_FP,
                  'rm': MctsPath.CAC_SMILES_NON_STEREO_RM_FP}

    target_inchi = None
    inchis = None
    smiles = None

    # def __init__(self, target_inchi: str = None, target_smiles: str = None, cac_types: [str] = None):
    #     print('loading cac db')
    #     self.target_inchi = target_inchi
    #     self._cac_types = ['bbm', 'sm', 'rm'] if cac_types is None else cac_types
    #     self._inchis = self.load_inchis()
    #     self._smiles = self.load_smiles()
    #     # self.remove_target(target_inchi, target_smiles)
    #     print('loaded cac db')

    @classmethod
    def init(cls, cac_types=None):
        if cac_types is None:
            cac_types = ['bbm', 'sm', 'rm']
        cls.inchis = cls.load_inchis(cac_types)
        cls.smiles = cls.load_smiles(cac_types)

    @classmethod
    def is_cac(cls, inchi: str = None, smiles: str = None) -> bool:
        if inchi == cls.target_inchi:
            return False
        if inchi is None and smiles is None:
            raise ValueError('cac: inchi and smiles both is None')
        res = False
        if inchi is not None and inchi in cls.inchis:
            res = True
        if smiles is not None and smiles in cls.smiles:
            res = True
        return res

    # ======================================

    @classmethod
    def _remove_inchi_stereo(cls, inchi: str):
        mol = AllChem.MolFromInchi(inchi)
        AllChem.RemoveStereochemistry(mol)
        return AllChem.MolToInchi(mol)

    @classmethod
    def _remove_smiles_stereo(cls, smiles: str):
        mol = AllChem.MolFromSmiles(smiles)
        AllChem.RemoveStereochemistry(mol)
        return AllChem.MolToSmiles(mol)

    # ==================================

    @classmethod
    def load_inchis(cls, cac_types):
        res = set()
        for cac_type in cac_types:
            res = res.union(cls._load_non_stereo_inchis_set(cls.inchi_fps[cac_type]))
        res = res.union(cls.load_stereo_inchis())
        res = res - {'InChI=1S/C14H13N5O4/c1-7-17-6-19(18-7)13-11-10(9(22-2)5-16-13)8(4-15-11)12(20)14(21)23-3/h4-6,15H,1-3H3',
                     'InChI=1S/C24H23N7O4/c1-15-27-14-31(28-15)22-20-19(18(35-2)13-26-22)17(12-25-20)21(32)24(34)30-10-8-29(9-11-30)23(33)16-6-4-3-5-7-16/h3-7,12-14,25H,8-11H2,1-2H3',
                     # 'InChI=1S/C11H11N5O/c1-7-14-6-16(15-7)11-10-8(3-4-12-10)9(17-2)5-13-11/h3-6,12H,1-2H3',
                     'InChI=1S/C21H20N4O2/c1-12(2)26-19-9-6-13(10-14(19)11-22)21-24-20(25-27-21)17-5-3-4-16-15(17)7-8-18(16)23/h3-6,9-10,12,18H,7-8,23H2,1-2H3',
                     'InChI=1S/H3NO/c1-2/h2H,1H2',
                     'InChI=1S/C16H24N2O/c17-10-8-15(14-5-1-4-11-18-14)9-12-19-16(13-15)6-2-3-7-16/h1,4-5,11H,2-3,6-10,12-13,17H2',
                     'InChI=1S/C25H36FN5O3S/c1-23(2,3)34-22(32)31-14-12-30(13-15-31)21-27-16-19(17-28-21)25(7,29-35(33)24(4,5)6)18-8-10-20(26)11-9-18/h8-11,16-17,29H,12-15H2,1-7H3',
                     }
        res.add('InChI=1S/C6H5IN4/c7-5-2-1-4-6(8)9-3-10-11(4)5/h1-3H,(H2,8,9,10)')
        res.add('InChI=1S/C15H23ClNO4P/c1-4-13(5-2)11-20-15(18)12(3)17-22(16,19)21-14-9-7-6-8-10-14/h6-10,12-13H,4-5,11H2,1-3H3,(H,17,19)')
        print(f'load cac inchis: {len(res)}')
        return res

    @classmethod
    def load_stereo_inchis(cls):
        res = set()
        bbm_df = pd.read_csv(MctsPath.CAC_INCHI_BBM_FP, sep='\t', encoding='utf-8')
        sm_df = pd.read_csv(MctsPath.CAC_INCHI_SM_FP, sep='\t', encoding='utf-8')
        rm_df = pd.read_csv(MctsPath.CAC_INCHI_RM_FP, sep='\t', encoding='utf-8')
        res = res.union(set(bbm_df['inchi']))
        res = res.union(set(sm_df['inchi']))
        res = res.union(set(rm_df['inchi']))
        return res

    @classmethod
    def load_smiles(cls, cac_types):
        res = set()
        for cac_type in cac_types:
            res = res.union(cls._load_smiles_set(cls.smiles_fps[cac_type]))
        print(f'load cac smiles: {len(res)}')
        return res

    # ========================

    @classmethod
    def _load_smiles_set(cls, fp):
        with open(fp, 'r', encoding='utf-8')as f:
            return set(f.read().split('\n'))

    @classmethod
    def _load_inchis_set(cls, fp):
        bbm_df = pd.read_csv(fp, sep='\t', encoding='utf-8')
        return set(bbm_df['inchi'])

    @classmethod
    def _load_non_stereo_inchis_set(cls, fp):
        with open(fp, 'r', encoding='utf-8')as f:
            return set(f.read().split('\n'))

    # ========================

    @classmethod
    def _convert_inchis_to_smiles(cls, inchis):
        for inchi in inchis:
            try:
                mol = AllChem.MolFromInchi(inchi)
            except:
                continue
            if mol is None:
                continue
            yield AllChem.MolToSmiles(mol)

    @classmethod
    def _save_smiles(cls, smiles_iter, fp):
        smileses = list(smiles_iter)
        with open(fp, 'w', encoding='utf-8')as f:
            f.write('\n'.join(smileses))

    @classmethod
    def get_smiles(cls):
        bbm_smiles_iter = cls._convert_inchis_to_smiles(cls._load_inchis_set(MctsPath.CAC_INCHI_BBM_FP))
        cls._save_smiles(bbm_smiles_iter, MctsPath.CAC_SMILES_BBM_FP)

        sm_smiles_iter = cls._convert_inchis_to_smiles(cls._load_inchis_set(MctsPath.CAC_INCHI_SM_FP))
        cls._save_smiles(sm_smiles_iter, MctsPath.CAC_SMILES_SM_FP)

        rm_smiles_iter = cls._convert_inchis_to_smiles(cls._load_inchis_set(MctsPath.CAC_INCHI_RM_FP))
        cls._save_smiles(rm_smiles_iter, MctsPath.CAC_SMILES_RM_FP)

        bbm_non_stereo_smiles_iter = cls._convert_inchis_to_smiles(cls._load_non_stereo_inchis_set(MctsPath.CAC_INCHI_NON_STEREO_BBM_FP))
        cls._save_smiles(bbm_non_stereo_smiles_iter, MctsPath.CAC_SMILES_NON_STEREO_BBM_FP)

        sm_non_stereo_smiles_iter = cls._convert_inchis_to_smiles(cls._load_non_stereo_inchis_set(MctsPath.CAC_INCHI_NON_STEREO_SM_FP))
        cls._save_smiles(sm_non_stereo_smiles_iter, MctsPath.CAC_SMILES_NON_STEREO_SM_FP)

        rm_non_stereo_smiles_iter = cls._convert_inchis_to_smiles(cls._load_non_stereo_inchis_set(MctsPath.CAC_INCHI_NON_STEREO_RM_FP))
        cls._save_smiles(rm_non_stereo_smiles_iter, MctsPath.CAC_SMILES_NON_STEREO_RM_FP)


if __name__ == "__main__":
    # Cac.get_smiles()
    # cac = DbCac('InChI=1S/C5H6Cl2O/c6-3-1-2-4(7)5(3)8/h3-4H,1-2H2')
    DbCac.init()
    # print(DbCac.is_cac('InChI=1S/C6H10O2/c1-4-5(2)8-6(3)7/h4-5H,1H2,2-3H3'))
    # print(DbCac.is_cac('InChI=1S/ClH/h1H'))
    print(DbCac.is_cac('InChI=1S/C15H15FN4O/c16-13-3-1-11(2-4-13)14(21)12-9-18-15(19-10-12)20-7-5-17-6-8-20/h1-4,9-10,17H,5-8H2'))
