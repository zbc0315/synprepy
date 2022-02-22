#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Time    : 2022/1/7 14:19
# @Author  : zhangbc0315@outlook.com
# @File    : rxn_templates_analyzer.py
# @Software: PyCharm

from matplotlib import pyplot as plt


from data_utils import RxnTemplates
from config import Config


class RxnTemplatesAnalyzer:

    def __init__(self, config: Config):
        self._rxn_centralized_templates = RxnTemplates(config.rxn_centralized_template_tsv_file_path)
        self._rxn_extended_templates = RxnTemplates(config.rxn_extended_template_tsv_file_path)

    @classmethod
    def _get_draw_data(cls, data: RxnTemplates):
        xs = []
        ys = []
        for x in range(100):
            xs.append(x)
            ys.append(data.calc_coverage(x))
        return xs, ys

    def draw(self):
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 4))
        cxs, cys = self._get_draw_data(self._rxn_centralized_templates)
        exs, eys = self._get_draw_data(self._rxn_extended_templates)
        ax1.plot(cxs, cys)
        ax2.plot(exs, eys)
        plt.show()


if __name__ == "__main__":
    RxnTemplatesAnalyzer(Config("../config.json")).draw()
