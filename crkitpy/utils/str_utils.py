# -*- coding: utf-8 -*-
# @Time     : 2020/4/28 11:09
# @Author   : zbc
# @Email    : zhangbc0315@outlook.com
# @FileName : str_utils.py
# @Software : PyCharm

from typing import Tuple


class StrUtils:

    @staticmethod
    def get_int_from_idx(string: str, idx: int) -> Tuple[int, int]:
        int_c_list = []
        for i in range(idx, len(string)):
            if string[i].isdigit():
                int_c_list.append(string[i])
            else:
                idx = i
                break
        if len(int_c_list) > 0:
            return int(''.join(int_c_list)), idx
        else:
            return None, idx
