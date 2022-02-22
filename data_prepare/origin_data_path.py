# -*- coding: utf-8 -*-
# @Time     : 2021/3/3 14:21
# @Author   : zbc
# @Email    : zhangbc0315@outlook.com
# @FileName : origin_data_path.py
# @Software : PyCharm

import os

from proj_path import PROJ_PATH

KNOWN_CAT_INCHIS = os.path.join(PROJ_PATH, 'Spectators/knownCatalystInChIs.txt')
KNOWN_SOL_INCHIS = os.path.join(PROJ_PATH, 'Spectators/knownSolventInChIs.txt')

ALL_RXN_FP = os.path.join(PROJ_PATH, 'chemdb/rxndb_public_clean_duplicate_rxn.tsv')
ALL_RXN_WITH_SPEC_FP = os.path.join(PROJ_PATH, 'chemdb/clean_duplicate_rxn_with_spec.tsv')


if __name__ == '__main__':
    pass
