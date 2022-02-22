# -*- coding: utf-8 -*-
# @Time     : 2021/3/30 15:26
# @Author   : zbc
# @Email    : zhangbc0315@outlook.com
# @FileName : rxn_mol_table.py
# @Software : PyCharm

from sqlalchemy import Column, Integer, ForeignKey, SmallInteger

from rmp_database.chem_db_orm.base_db import Base
from rmp_database.chem_db_orm.rxn_table import RxnTable
from rmp_database.chem_db_orm.mol_table import MolTable


class RxnMolTable(Base):

    __tablename__ = 'rxn_mol'

    rmid = Column(Integer, primary_key=True, autoincrement=True)
    # 分子在反应中的角色：0-反应物，1-催化剂，2-溶剂，3-产物, 4-其他
    mol_role = Column(SmallInteger)
    rid = Column(Integer, ForeignKey("rxn.rid"))
    mid = Column(Integer, ForeignKey("mol.mid"))


if __name__ == '__main__':
    pass
