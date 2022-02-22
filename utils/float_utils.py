#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Time    : 2022/2/12 14:32
# @Author  : zhangbc0315@outlook.com
# @File    : float_utils.py
# @Software: PyCharm


class FloatUtils:

    @classmethod
    def to_str(cls, f: float, p: int) -> str:
        s = str(f)
        if '.' not in s:
            return s
        else:
            integer, decimal = s.split('.')
            if len(decimal) <= p:
                return s
            else:
                decimal = decimal[:p]
                return f"{integer}.{decimal}"


if __name__ == "__main__":
    pass
