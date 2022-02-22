# -*- coding: utf-8 -*-
# @Time     : 2021/4/26 16:22
# @Author   : zbc
# @Email    : zhangbc0315@outlook.com
# @FileName : known_spectators_loader.py
# @Software : PyCharm

from typing import Iterator

from rdkit.Chem import AllChem

from config import SSPath
from log_utils import Logger


class KnownSpectatorsLoader:

    @classmethod
    def _load_spectators_inchis(cls, fp: str) -> Iterator[str]:
        """

        :param fp:
        :return:
        """
        with open(fp, 'r', encoding="utf-8")as f:
            for line in f.readlines():
                line = line.strip()
                if line.startswith('#') or len(line) == 0:
                    continue
                yield line.split('\t')[0]

    @classmethod
    def load_cat_inchis(cls) -> Iterator[str]:
        for inchi in cls._load_spectators_inchis(SSPath.KNOWN_CAT_INCHIS):
            yield inchi

    @classmethod
    def load_cat_smiles(cls) -> Iterator[str]:
        for inchi in cls.load_cat_inchis():
            try:
                yield AllChem.MolToSmiles(AllChem.MolFromInchi(inchi))
            except Exception as e:
                Logger.warn(f"Wrong Inchi: {inchi} with error {e}")

    @classmethod
    def load_sol_inchis(cls) -> Iterator[str]:
        for inchi in cls._load_spectators_inchis(SSPath.KNOWN_SOL_INCHIS):
            yield inchi

    @classmethod
    def load_sol_smiles(cls):
        for inchi in cls.load_sol_inchis():
            yield AllChem.MolToSmiles(AllChem.MolFromInchi(inchi))


if __name__ == '__main__':
    pass
