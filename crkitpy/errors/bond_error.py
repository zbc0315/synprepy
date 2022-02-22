# -*- coding: utf-8 -*-
# @Time     : 2020/6/11 15:40
# @Author   : zbc
# @Email    : zhangbc0315@outlook.com
# @FileName : bond_error.py
# @Software : PyCharm


class BondError(Exception):

    def __init__(self, *args):
        self.args = args


class BondTypeError(BondError):

    def __init__(self, *args):
        self.args = args


class BondValenceError(BondError):

    def __init__(self, *args):
        self.args = args
