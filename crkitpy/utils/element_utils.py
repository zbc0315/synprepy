# -*- coding: utf-8 -*-
# @Time     : 2020/4/23 14:59
# @Author   : zbc
# @Email    : zhangbc0315@outlook.com
# @FileName : element_utils.py
# @Software : PyCharm

symbol_to_color = {'C': 'black', 'N': 'blue', 'O': 'red',
                   'P': 'brown', 'S': 'brown', 'Cl': 'green'}


class ElementUtils:

    @classmethod
    def get_color_by_symbol(cls, symbol: str) -> str:
        return symbol_to_color[symbol]

