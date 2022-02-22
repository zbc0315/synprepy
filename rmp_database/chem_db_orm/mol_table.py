# -*- coding: utf-8 -*-
# @Time     : 2021/3/30 15:25
# @Author   : zbc
# @Email    : zhangbc0315@outlook.com
# @FileName : mol_table.py
# @Software : PyCharm

from sqlalchemy import Column, Integer, Text
from sqlalchemy.orm import relationship

from rmp_database.chem_db_orm.base_db import Base


class MolTable(Base):

    __tablename__ = 'mol'

    mid = Column(Integer, primary_key=True, autoincrement=True)
    name_ch = Column(Text)
    name_en = Column(Text)
    smiles = Column(Text)
    inchi = Column(Text)
    rxn_mol = relationship("RxnMolTable", backref="mol")


if __name__ == '__main__':
    pass
