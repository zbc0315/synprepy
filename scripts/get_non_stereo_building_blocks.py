# -*- coding: utf-8 -*-
# @Time     : 2021/5/8 14:40
# @Author   : zbc
# @Email    : zhangbc0315@outlook.com
# @FileName : get_non_stereo_building_blocks.py
# @Software : PyCharm

import pandas as pd
from rdkit.Chem import AllChem

from config.path import BuildingBlocksPath


class GetNonStereoBuildingBlocks:

    @classmethod
    def _remove_stereo(cls, inchi: str):
        mol = AllChem.MolFromInchi(inchi)
        AllChem.RemoveStereochemistry(mol)
        return AllChem.MolToInchi(mol)

    @classmethod
    def _get_non_stereo_inchi(cls, inchi: str):
        try:
            new_inchi = cls._remove_stereo(inchi)
        except Exception as e:
            new_inchi = None
        if new_inchi is None or new_inchi == inchi:
            return inchi, False
        return new_inchi, True

    @classmethod
    def _write_non_stereo_inchis(cls, inchis_fp: str, non_stereo_inchis_fp: str):
        non_stereo_inchis = []
        df = pd.read_csv(inchis_fp, sep='\t', encoding='utf-8')
        for n, inchi in enumerate(df["inchi"]):
            if n % 1000 == 0:
                print(f"{n} - {inchis_fp}")
            new_inchi, removed_stereo = cls._get_non_stereo_inchi(inchi)
            if not removed_stereo:
                continue
            if new_inchi not in non_stereo_inchis:
                non_stereo_inchis.append(new_inchi)
        with open(non_stereo_inchis_fp, 'w', encoding='utf-8')as f:
            f.write('\n'.join(non_stereo_inchis))

    @classmethod
    def process(cls):
        cls._write_non_stereo_inchis(BuildingBlocksPath.BBM_FP, BuildingBlocksPath.BBM_NON_STEREO_FP)
        cls._write_non_stereo_inchis(BuildingBlocksPath.SM_FP, BuildingBlocksPath.SM_NON_STEREO_FP)
        cls._write_non_stereo_inchis(BuildingBlocksPath.RM_FP, BuildingBlocksPath.RM_NON_STEREO_FP)


if __name__ == '__main__':
    GetNonStereoBuildingBlocks.process()
