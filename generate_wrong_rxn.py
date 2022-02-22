#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Time    : 2021/11/12 16:22
# @Author  : zhangbc0315@outlook.com
# @File    : generate_wrong_rxn.py
# @Software: PyCharm

import sys
from rxn_determiner.prepare_wrong_rxn import PrepareWrongRxn


if __name__ == "__main__":
    pid = int(sys.argv[1])
    PrepareWrongRxn.process(pid)
