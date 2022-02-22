# -*- coding: utf-8 -*-
# @Time     : 2021/3/24 14:52
# @Author   : zbc
# @Email    : zhangbc0315@outlook.com
# @FileName : logger.py
# @Software : PyCharm

import logging

from logging.handlers import TimedRotatingFileHandler

from mcts_lite.mcts_path import MctsPath


class Logging:

    log = None

    @classmethod
    def init(cls, process_id: int, level=logging.INFO):

        cls.log = logging.getLogger()
        cls.log.setLevel(level)
        filehandler = logging.handlers.TimedRotatingFileHandler(MctsPath.LOG_FP % process_id, when='d', interval=2,
                                                                backupCount=7, encoding='utf-8')
        formatter = logging.Formatter('%(levelname)s: %(asctime)s %(filename)s %(message)s')
        filehandler.suffix = "%Y-%m-%d_%H-%M-%S.log"
        filehandler.setFormatter(formatter)
        cls.log.addHandler(filehandler)
        cls.log.info('init log')

    @classmethod
    def get_log(cls):
        return cls.log


if __name__ == '__main__':
    logging.debug('debug')
