#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Time    : 2021/9/15 10:24
# @Author  : zhangbc0315@outlook.com
# @File    : mcts_evaluate_process.py
# @Software: PyCharm
import sys
from multiprocessing import Pool

from mcts_evaluate.mcts_evaluate_process import MctsEvaluate


def eva_process(pid):
    MctsEvaluate.process(pid)
    # try:
    #     MctsEvaluate.process(pid)
    # except Exception as e:
    #     print('='*100)
    #     print(f'[ERROR] Process {pid}')
    #     print(e)
    #     print('=' * 100)


def multi_process():
    num_pool = 20
    p = Pool(num_pool)
    for i in range(num_pool):
        p.apply_async(eva_process, args=(i, ))
        print(f'Begin Process: {i}')
    p.close()
    p.join()


if __name__ == "__main__":
    # MctsEvaluate.process(0)
    multi_process()
