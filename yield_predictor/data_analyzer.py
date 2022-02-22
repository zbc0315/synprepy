#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Time    : 2022/2/12 13:37
# @Author  : zhangbc0315@outlook.com
# @File    : data_analyzer.py
# @Software: PyCharm

import logging

import matplotlib.pyplot as plt

from utils import FloatUtils
from config import Config
from data_utils import MolData, RxnData, RxnDataType, SpectatorData, SpectatorType


class DataAnalyzer:

    def __init__(self, gap: float = 5):
        self._gap = gap
        self._config = Config("../config.json")
        self._mol_data = MolData(self._config.mol_data_tsv_file_path)
        self._rxn_data = RxnData(self._config.rxn_data_tsv_file_path,
                                 self._mol_data,
                                 RxnDataType.TRAIN_TEST,
                                 self._config.evaluate_year)
        self._rxn_data.init_yield()

    @classmethod
    def duplicate_rare_data(cls, rxn_df):
        start_end_count_rxn_codes = cls.count_rxn_code_yield(rxn_df)
        counts = [r[2] for r in start_end_count_rxn_codes]
        max_count = max(counts)
        res = []
        for start, end, count, rxn_codes in start_end_count_rxn_codes:
            ratio = int(max_count / count)
            if ratio <= 1:
                ratio = 1
            logging.info(f"Duplicate Rare Rxn Data yield({start}-{end}) from {count} to {count*ratio}")
            res.extend(rxn_codes * ratio)
        logging.info(f"Duplicate Rxn Data, get {len(res)} from {sum(counts)}")
        return res

    @classmethod
    def count_rxn_code_yield(cls, rxn_df):
        res = []
        start = 0
        while start < 100:
            end = min(start + 5, 100)
            rxn_codes = list(rxn_df.query(f"{start}<yield_product<={end}").rxn_code)
            res.append((start, end, len(rxn_codes), rxn_codes))
            start = end
        return res

    def count_yield(self):
        start = 0
        x_labels = []
        counts = []
        # start_end_count_rxn_codes = self.duplicate_rare_data(rxn_df=)
        while start < 100:
            end = min(start + self._gap, 100)
            x_labels.append(f"{FloatUtils.to_str(start, 1)}-{FloatUtils.to_str(end, 1)}")
            s = '<=' if start == 0 else '<'
            counts.append(len(self._rxn_data.rxn_df.query(f"{start}<yield_product<={end}")) / 1000)
            start += self._gap
        # plt.rcParams['figure.figsize'] = (6.0, 6.5)
        plt.rcParams['figure.dpi'] = 300
        plt.title("yield count")
        plt.xlabel("yield / %")
        plt.ylabel("count / 1000")
        plt.bar(range(len(counts)), counts, tick_label=x_labels)
        plt.xticks(rotation=270)
        for n, count in enumerate(counts):
            count_str = str(int(count))
            plt.text(n - 0.5, count - 1.5 * len(count_str), count_str, color="white", rotation=270)
        plt.tight_layout()
        plt.show()

    # def count_yield(self):
    #     start = 0
    #     x_labels = []
    #     counts = []
    #     while start < 100:
    #         end = min(start + self._gap, 100)
    #         x_labels.append(f"{FloatUtils.to_str(start, 1)}-{FloatUtils.to_str(end, 1)}")
    #         counts.append(len(self._rxn_data.rxn_df.query(f"{start}<yield_product<={end}"))/1000)
    #         start += self._gap
    #     # plt.rcParams['figure.figsize'] = (6.0, 6.5)
    #     plt.rcParams['figure.dpi'] = 300
    #     plt.title("yield count")
    #     plt.xlabel("yield / %")
    #     plt.ylabel("count / 1000")
    #     plt.bar(range(len(counts)), counts, tick_label=x_labels)
    #     plt.xticks(rotation=270)
    #     for n, count in enumerate(counts):
    #         count_str = str(int(count))
    #         plt.text(n-0.5, count-1.5*len(count_str), count_str, color="white", rotation=270)
    #     plt.tight_layout()
    #     plt.show()


if __name__ == "__main__":
    DataAnalyzer().count_yield()
