# -*- coding: utf-8 -*-
# @Time     : 2020/6/4 16:48
# @Author   : zbc
# @Email    : zhangbc0315@outlook.com
# @FileName : rxn_error.py
# @Software : PyCharm


class RxnError(Exception):

    def __init__(self, *args):
        self.args = args


class OneStepRxnError(RxnError):

    def __init__(self, *args):
        self.args = args


class ProductFullyMapError(RxnError):

    def __init__(self, *args):
        self.args = args


class ReactantUniqueMapNumError(RxnError):

    def __init__(self, *args):
        self.args = args


class ReactantDuplicateError(RxnError):

    def __init__(self, *args):
        self.args = args


class MoleculeNumError(RxnError):

    def __init__(self, *args):
        self.args = args


class TemplateGroupDuplicateError(RxnError):

    def __init__(self, *args):
        self.args = args
