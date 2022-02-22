# -*- coding: utf-8 -*-
# @Time     : 2020/6/19 10:07
# @Author   : zbc
# @Email    : zhangbc0315@outlook.com
# @FileName : smiles_error.py
# @Software : PyCharm


class SmilesError(Exception):

    def __init__(self, *args):
        self.args = args


class RdkitSmilesError(SmilesError):

    def __init__(self, *args):
        self.args = args


class SmilesWriterError(SmilesError):

    def __init__(self, *args):
        self.args = args
