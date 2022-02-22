#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Time    : 2021/9/13 11:20
# @Author  : zhangbc0315@outlook.com
# @File    : mcts_evaluate_process.py
# @Software: PyCharm
import datetime
import json
import logging
import os
from datetime import datetime
import time

import pandas as pd

from mcts_lite.mcts import Mcts
from mcts_lite.mcts_path import MctsPath
from mcts_lite.node_utils import NodeUtils
from log_utils.logger import Logging


class MctsEvaluate:

    @classmethod
    def get_products(cls):
        df = pd.read_csv(MctsPath.TARGET_PRODUCTS_FP, sep='\t', encoding='utf-8')
        products = set()
        for _, row in df.iterrows():
            if row['inchi'] in products:
                continue
            products.add(row['inchi'])
            yield row['rid'], row['inchi']

    @classmethod
    def get_unsolved_products_in10(cls):
        for fn in os.listdir(MctsPath.EVA10_DP):
            fp = os.path.join(MctsPath.EVA10_DP, fn)
            df = pd.read_csv(fp, sep='\t', encoding='utf-8')
            for _, row in df.iterrows():
                if not row.succeed:
                    yield row.rid, row.inchi

    @classmethod
    def _mcts(cls, inchi):
        rxn_path = None
        used_seconds = -1
        begin_time = datetime.now()
        for node, is_right in Mcts.search(inchi):
            if not is_right:
                continue
            rxn_path, _, _ = NodeUtils.get_rxn_plan(node)
            used_seconds = (datetime.now() - begin_time).seconds
            break
        return rxn_path, used_seconds

    @classmethod
    def write_line(cls, rid, succeed, used_seconds, rxn_path, inchi, result_fp):
        if rxn_path is None:
            rxn_path = ''
        else:
            rxn_path = json.dumps(rxn_path)
        with open(result_fp, 'a', encoding='utf-8') as f:
            f.write(f'\n{rid}'
                    f'\t{succeed}'
                    f'\t{used_seconds}'
                    f'\t{rxn_path}'
                    f'\t{inchi}')

    @classmethod
    def init_file(cls, result_fp):
        with open(result_fp, 'w', encoding='utf-8') as f:
            f.write('rid\tsucceed\tused_seconds\trxn_path\tinchi')

    @classmethod
    def process(cls, pid):
        time.sleep(pid)
        Logging.init(pid+300, logging.WARN)
        Logging.log.info(f'PROCESS ID: {pid}')
        Logging.log.warning(f'Test Warning Evaluate: {pid}')
        result_fp = MctsPath.EVALUATE30_RESULT_FP % pid
        if os.path.exists(result_fp):
            saved_df = pd.read_csv(result_fp, sep='\t', encoding='utf8')
        else:
            cls.init_file(result_fp)
            saved_df = None
        succeed_num = 0
        failed_num = 0
        Mcts.init(25, 1800, 1)
        # for n, (rid, inchi) in enumerate(cls.get_products()):
        for n, (rid, inchi) in enumerate(cls.get_unsolved_products_in10()):
            if n % 20 != pid:
                continue
            if saved_df is not None and rid in list(saved_df.rid):
                Logging.log.info(f'parsed {n} - {rid} - {inchi}')
                continue
            Logging.log.info(inchi)
            try:
                rxn_path, used_seconds = cls._mcts(inchi)
            except Exception as e:
                Logging.log.info('#'*100)
                Logging.log.info(f'error inchi: {inchi}')
                Logging.log.info(e)
                Logging.log.info('#'*100)
                continue
            cls.write_line(rid, rxn_path is not None, used_seconds, rxn_path, inchi, result_fp)
            # res['rid'].append(rid)
            # res['succeed'].append(rxn_path is not None)
            # res['used_seconds'].append(used_seconds)
            # res['inchi'].append(inchi)
            # res['rxn_path'].append(json.dumps(rxn_path))
            if rxn_path is None:
                failed_num += 1
            else:
                succeed_num += 1
            Logging.log.info(f'succeed: {succeed_num}, failed: {failed_num}')
            print(f'succeed: {succeed_num}, failed: {failed_num}')
        # df = pd.DataFrame(res)
        # df.to_csv(MctsPath. EVALUATE_RESULT_FP, sep='\t', encoding='utf-8', index=False)


if __name__ == "__main__":
    MctsEvaluate.process(0)
