# -*- coding: utf-8 -*-
# @Time     : 2020/6/5 14:54
# @Author   : zbc
# @Email    : zhangbc0315@outlook.com
# @FileName : atom_error.py
# @Software : PyCharm


class AtomError(Exception):

    def __init__(self, *args):
        self.args = args


class AtomValenceError(AtomError):

    def __init__(self, *args):
        self.args = args
