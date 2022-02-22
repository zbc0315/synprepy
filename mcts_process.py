#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Time    : 2021/9/16 12:48
# @Author  : zhangbc0315@outlook.com
# @File    : mcts_process.py
# @Software: PyCharm

import sys

from synpre_backend.mcts_process import ProcessRunner
from log_utils.logger import Logging


def run_process(pid):
    Logging.init(pid)
    Logging.log.info(f'Process: {pid}')
    ProcessRunner().process(pid)


def run(pid):

    while True:
        run_process(pid)


if __name__ == "__main__":
    run_process(int(sys.argv[1]))

