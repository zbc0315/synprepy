#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Time    : 2021/8/12 22:20
# @Author  : zhangbc0315@outlook.com
# @File    : move_rxn_data_to_mysql.py
# @Software: PyCharm

from rmp_database.old_chem_db import OldChemDb

from rmp_database.chem_db_orm import get_session
from rmp_database.chem_db_orm import MolTable, RxnTable, RxnMolTable


class MoveRxnDataToMysql:

    @classmethod
    def get_rxn_data_before_2016(cls):
        odb = OldChemDb()
        for rxn_data in odb.get_data_iter("public.clean_duplicate_rxn",
                                          ['rid', 'rxn_smi', 'catalysts_sids', 'solvents_sids'],
                                          None):
            pass

    @classmethod
    def add_rxn_mol(cls):
        rxn = RxnTable()

    @classmethod
    def process_before_2016(cls):
        pass


if __name__ == "__main__":
    pass
