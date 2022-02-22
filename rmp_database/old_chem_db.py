# -*- coding: utf-8 -*-
# @Time     : 2021/3/3 13:59
# @Author   : zbc
# @Email    : zhangbc0315@outlook.com
# @FileName : old_chem_db.py
# @Software : PyCharm

# class OldChemDb:
#
#     def __init__(self):
#         pass

from typing import Union, Iterator, Dict

import psycopg2
from psycopg2 import extensions


class OldChemDb:

    def __init__(self):
        self.conn: extensions.connection = self._create_connect()
        self.cur: extensions.cursor = self.conn.cursor()

    @staticmethod
    def _create_connect() -> extensions.connection:
        conn = psycopg2.connect(database='rxndb',
                                user='rxnpredictor',
                                password='password',
                                host='127.0.0.1',
                                port='5444')
        return conn

    def commit(self):
        self.conn.commit()

    @staticmethod
    def _get_value(v):
        if isinstance(v, str):
            return f"'{v}'"
        return v

    @classmethod
    def _get_update_values(cls, cols: [str], values: []) -> str:
        res_list = []
        for col, value in zip(cols, values):
            res_list.append(f"{col} = {cls._get_value(value)}")
        return ",".join(res_list)

    def update(self, table_name:str, cols: [str], values: [], condition: str, commit: bool = True):
        sql = f"UPDATE {table_name} SET {self._get_update_values(cols, values)} WHERE {condition}"
        self.cur.execute(sql)
        if commit:
            self.commit()

    def insert_one_data(self, table_name: str, cols: [str], values: [], commit: bool = True):
        cols_str = ",".join(cols)
        vals_str = ",".join([str(self._get_value(v)) for v in values])
        sql = f"INSERT INTO {table_name} ({cols_str}) VALUES ({vals_str})"
        self.cur.execute(sql)
        if commit:
            self.commit()

    def search_one_data(self, table_name: str, cols: [str], conditions: Union[str, None]) -> Union[Dict, None]:
        self.cur.execute(self._select_sql(table_name, cols, conditions, None, False))
        try:
            res = self.cur.fetchone()
            res_dict = {}
            for col, value in zip(cols, res):
                res_dict[col] = value
            return res_dict
        except Exception as e:
            return None

    def get_data_iter(self, table_name: str, cols: [str], conditions: Union[str, None], order_by: str = None, desc: bool = False) -> Iterator[Dict]:
        self.cur.execute(self._select_sql(table_name, cols, conditions, order_by, desc))
        while True:
            try:
                res = self.cur.fetchone()
                res_dict = {}
                for col, value in zip(cols, res):
                    res_dict[col] = value
                yield res_dict
            except Exception as e:
                break

    @staticmethod
    def _select_sql(table_name: str, cols: [str], condition: str, order_by: str, desc: bool) -> str:
        condition_sql = f" WHERE {condition}" if condition is not None else ""
        desc_str = " DESC" if desc else " ASC"
        order_by_sql = f" ORDER BY {order_by}{desc_str}" if order_by is not None else ""
        return f"SELECT {','.join(cols)} FROM {table_name}{condition_sql}{order_by_sql}"
        if condition is not None:
            return f"SELECT {','.join(cols)} " \
                   f"FROM {table_name} " \
                   f"WHERE {condition}"
        else:
            return f"SELECT {','.join(cols)} " \
                   f"FROM {table_name} "


if __name__ == '__main__':
    pass
