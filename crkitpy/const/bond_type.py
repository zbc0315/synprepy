# -*- coding: utf-8 -*-
# @Time     : 2020/4/20 16:32
# @Author   : zbc
# @Email    : zhangbc0315@outlook.com
# @FileName : bond_type.py
# @Software : PyCharm

SINGLE = 'SINGLE'
DOUBLE = 'DOUBLE'
TRIPLE = 'TRIPLE'
AROMATIC = 'AROMATIC'

BOND_WEIGHTS = {SINGLE: 1,
                DOUBLE: 2,
                TRIPLE: 3,
                AROMATIC: 1.5}

BOND_WEIGHTS_TO_TYPE = {
    1: SINGLE,
    2: DOUBLE,
    3: TRIPLE,
    1.5: AROMATIC
}

BOND_WEIGHT_STR_TO_SYMBOLS = {'1': SINGLE,
                              '2': DOUBLE,
                              '3': TRIPLE,
                              '1.5': AROMATIC}

BOND_SYMBOLS = {'-': SINGLE,
                '=': DOUBLE,
                '#': TRIPLE,
                ':': AROMATIC}

BOND_TO_SYMBOLS = {SINGLE: '-',
                   DOUBLE: '=',
                   TRIPLE: '#',
                   AROMATIC: ':'}
