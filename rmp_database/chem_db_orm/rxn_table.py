# -*- coding: utf-8 -*-
# @Time     : 2021/3/30 15:25
# @Author   : zbc
# @Email    : zhangbc0315@outlook.com
# @FileName : rxn_table.py
# @Software : PyCharm

from sqlalchemy import Column, Integer, String, Text, Date
from sqlalchemy.orm import relationship

from rmp_database.chem_db_orm.base_db import Base


class RxnTable(Base):

    __tablename__ = "rxn"

    rid = Column(Integer, primary_key=True, autoincrement=True)
    rxn_mol = relationship("RxnMolTable", backref='rxn')
    # code: 反应物id>>产物id,催化剂id,溶剂id
    rxn_code = Column(String(200))
    other_code = Column(String(100))
    rxn_smiles = Column(Text)
    source = Column(String(100))
    title = Column(Text)
    paragraph = Column(Text)
    year = Column(Integer)
    month = Column(Integer)
    day = Column(Integer)


if __name__ == '__main__':
    pass
