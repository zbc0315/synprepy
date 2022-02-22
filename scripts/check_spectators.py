# -*- coding: utf-8 -*-
# @Time     : 2021/4/26 16:21
# @Author   : zbc
# @Email    : zhangbc0315@outlook.com
# @FileName : check_spectators.py
# @Software : PyCharm

from log_utils import Logger
from data_utils.known_spectators_loader import KnownSpectatorsLoader
from database.old_chem_db import OldChemDb


class CheckSpectators:

    _known_cats_inchis = []
    _known_cats_smiles = []

    _known_sols_inchis = []
    _known_sols_smiles = []

    _unsaved_known_cats_num = 0
    _unsaved_known_sols_num = 0

    _odb = OldChemDb()

    @classmethod
    def load_known_spectators(cls) -> None:
        cls._known_cats_inchis = list(KnownSpectatorsLoader.load_cat_inchis())
        cls._known_cats_smiles = list(KnownSpectatorsLoader.load_cat_smiles())
        Logger.info(f"load known catalysts: {len(cls._known_cats_smiles)}")
        cls._known_sols_inchis = list(KnownSpectatorsLoader.load_sol_inchis())
        cls._known_sols_smiles = list(KnownSpectatorsLoader.load_sol_smiles())
        Logger.info(f"load known solvents: {len(cls._known_sols_smiles)}")

    @classmethod
    def _search_spectators(cls, smiles: str, inchi: str, role: str):
        spectator = cls._odb.search_one_data("public.spectators_mols",
                                             ["been_catalyst_times", "been_solvent_times"],
                                             f"inchi='{inchi}'")
        if spectator is None:
            if role == "cat":
                cls._unsaved_known_cats_num += 1
                num = cls._unsaved_known_cats_num
            else:
                cls._unsaved_known_sols_num += 1
                num = cls._unsaved_known_sols_num
            print(f"NO SAVED {role} {num} {smiles}")

    @classmethod
    def _process_cat_or_sol(cls, inchi_list: [str], smiles_list: [str], role: str):
        for smiles, inchi in zip(smiles_list, inchi_list):
            cls._search_spectators(smiles, inchi, role)

    @classmethod
    def process(cls):
        cls.load_known_spectators()
        cls._process_cat_or_sol(cls._known_cats_inchis, cls._known_cats_smiles, "cat")
        cls._process_cat_or_sol(cls._known_sols_inchis, cls._known_sols_smiles, "sol")


if __name__ == '__main__':
    CheckSpectators.process()
