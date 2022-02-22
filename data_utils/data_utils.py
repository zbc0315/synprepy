# -*- coding: utf-8 -*-
# @Time     : 2021/4/16 14:47
# @Author   : zbc
# @Email    : zhangbc0315@outlook.com
# @FileName : data_utils.py
# @Software : PyCharm

from typing import Tuple, List

import pandas as pd

from rmp_database.old_chem_db import OldChemDb
from data_utils.tid_utils import TidUtils


class DataUtils:

    @classmethod
    def iter_db_data(cls, data_tn: str, pid: int, processed_idx: int) -> Tuple[int, str, List[float]]:
        """

        :return: smiles, encoded tid
        """
        df = pd.read_csv(data_tn, sep='\t', encoding='utf-8')
        for idx, row in df.iterrows():
            if idx < processed_idx:
                continue
            # if idx % 10 != pid:
            #     continue
            yield idx, row.product_smi_with_map, TidUtils.encode_tids([row.tid])
        # db = OldChemDb()
        # for mol_data in db.get_data_iter(data_tn, ['product_smi_with_map', 'tid'], None):
        #     try:
        #         yield mol_data['product_smi_with_map'], TidUtils.encode_tids([mol_data['tid']])
        #     except KeyError as e:
        #         print(f"error tid: {mol_data['tid']}")
        #         yield None, None

    @classmethod
    def get_all_spectators_after_2016(cls) -> ({int}, {int}):
        """ 从数据库中获得

        :return:
        """
        pass

    @classmethod
    def get_all_spectators_before_2016(cls) -> ({int}, {int}):
        """

        :return:
        """
        pass

    @classmethod
    def get_all_spectators(cls) -> ({int}, {int}):
        """ 从已有数据中返回所有催化剂、溶剂的id

        :return: 催化剂id列表，溶剂id列表
        """
        pass


if __name__ == '__main__':
    pass
