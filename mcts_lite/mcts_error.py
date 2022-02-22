#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Time    : 2021/9/14 11:04
# @Author  : zhangbc0315@outlook.com
# @File    : mcts_error.py
# @Software: PyCharm


class ExpectedError(Exception):

    def __init__(self, *args):
        self.args = args


class UnexpectedError(Exception):

    def __init__(self, *args):
        self.args = args


if __name__ == "__main__":
    pass
