# -*- coding: utf-8 -*-
# @Time     : 2021/3/30 16:05
# @Author   : zbc
# @Email    : zhangbc0315@outlook.com
# @FileName : create_all.py
# @Software : PyCharm


from rmp_database.chem_db_orm.base_db import engine, Base
from rmp_database.chem_db_orm.mol_table import MolTable
from rmp_database.chem_db_orm.rxn_table import RxnTable
from rmp_database.chem_db_orm.rxn_mol_table import RxnMolTable


if __name__ == '__main__':
    Base.metadata.create_all(engine)
