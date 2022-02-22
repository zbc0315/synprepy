#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Time    : 2021/9/27 16:40
# @Author  : zhangbc0315@outlook.com
# @File    : link_spectators_smiles_to_sid.py
# @Software: PyCharm

import pandas as pd
from rdkit.Chem import AllChem

from config import SSPath
from rmp_database.old_chem_db import OldChemDb

from rdkit import RDLogger

RDLogger.DisableLog('rdApp.*')


class LinkSpectatorsSmilesToSid:

    odb = OldChemDb()

    @classmethod
    def get_sid_from_smiles(cls, smiles: str) -> int:
        """
        Args:
            smiles
        Return:
            如果查询到sid，则返回正数
            查询不到sid，返回-1
        """
        try:
            inchi = AllChem.MolToInchi(AllChem.MolFromSmiles(smiles))
        except:
            return -1
        s_data = cls.odb.search_one_data('public.spectators_mols', ['sid'], f"inchi='{inchi}'")
        return s_data['sid'] if s_data is not None else -1

    @classmethod
    def _process_fp(cls, source_fp: str, target_fp):
        df = pd.read_csv(source_fp, sep='\t', encoding='utf-8')
        df['sid'] = 0
        for i, row in df.iterrows():
            sid = cls.get_sid_from_smiles(row['smiles'])
            df.loc[i, 'sid'] = sid
        df.to_csv(target_fp, sep='\t', encoding='utf-8', index=False)

    @classmethod
    def process(cls):
        cls._process_fp(SSPath.CATS_COUNT_FP, SSPath.CATS_COUNT_WITH_SID_FP)
        cls._process_fp(SSPath.SOLS_COUNT_FP, SSPath.SOLS_COUNT_WITH_SID_FP)


if __name__ == "__main__":
    # print(LinkSpectatorsSmilesToSid.get_sid_from_smiles('[Ag-]'))
    LinkSpectatorsSmilesToSid.process()
