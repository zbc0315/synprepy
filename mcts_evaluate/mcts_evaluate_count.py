#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Time    : 2021/9/15 11:32
# @Author  : zhangbc0315@outlook.com
# @File    : mcts_evaluate_count.py
# @Software: PyCharm
import os

import pandas as pd
from matplotlib import pyplot as plt

from mcts_lite.mcts_path import MctsPath


class MctsEvaluateCount:

    @classmethod
    def count_error_inchis(cls):
        n = 0
        for fn in os.listdir(MctsPath.LOG_DP):
            fp = os.path.join(MctsPath.LOG_DP, fn)
            with open(fp, 'r', encoding='utf-8')as f:
                for line in f.readlines():
                    if 'Can not expand from root' in line:
                        n += 1
        print(f'error: {n}')
        return n

    @classmethod
    def get_eva_df(cls):
        res_df = None
        # fps = os.path.join('C:\\Users\\zhang\\Documents\\Data\\RxnMolPredictor\\Mcts\\Evaluate_', 'evaluate_result_%d.tsv')
        fps = MctsPath.EVALUATE_RESULT_FP
        for i in range(20):
            df_fp = fps % i
            if not os.path.exists(df_fp):
                print(f'No Exists: {df_fp}')
                continue
            df = pd.read_csv(df_fp, sep='\t', encoding='utf-8')
            res_df = df if res_df is None else res_df.append(df)
        return res_df

    @classmethod
    def _draw_time_used(cls, df):
        max_time = df['used_seconds'].max()
        max_time = 3600
        times = list(range(5, max_time+1))
        counts = []
        max_num = len(df)
        df = df.query('succeed')
        # max_num = len(df)
        num = 0
        for t in times:
            num = len(df.query(f'used_seconds<={t}'))
            counts.append(num/max_num)
        plt.plot(times, counts)
        plt.xlabel('Time/s')
        plt.ylabel('Solved Ratio')
        plt.show()
        print(f'max_num: {max_num}, right_num: {num}, ratio: {num/max_num}')

    @classmethod
    def process(cls):
        error_num = cls.count_error_inchis()
        error_num = 0
        df = cls.get_eva_df()
        count_all = len(df) - error_num
        right_df = df.query(f'succeed')
        count_right = len(right_df)
        mt = right_df.used_seconds.mean()
        it = right_df.used_seconds.median()
        print(f'{count_right/count_all} r({count_right}) e({count_all - count_right}) a({count_all}) mt({mt}) it({it})')
        cls._draw_time_used(df)


if __name__ == "__main__":
    MctsEvaluateCount.process()
