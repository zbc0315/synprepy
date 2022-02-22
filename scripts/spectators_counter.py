# -*- coding: utf-8 -*-
# @Time     : 2021/6/15 9:56
# @Author   : zbc
# @Email    : zhangbc0315@outlook.com
# @FileName : spectators_counter.py
# @Software : PyCharm

import pandas as pd

from config import SSPath


class SpectatorsCounter:

    _spec_count_fps = {"cat": SSPath.CATS_COUNT_FP,
                       "sol": SSPath.SOLS_COUNT_FP}
    _spec_times_count_fps = {"cat": SSPath.CATS_TIMES_COUNT_FP,
                             "sol": SSPath.SOLS_TIMES_COUNT_FP}

    @classmethod
    def _write_times_count(cls, sort_times_count: pd.Series, fp: str):
        with open(fp, 'w', encoding="utf-8")as f:
            f.write("times\tcount\n")
            for (time, count) in sort_times_count.items():
                f.write(f"{time}\t{count}\n")

    @classmethod
    def _count_spec(cls, spec_type: str):
        spec_fp = cls._spec_count_fps[spec_type]
        spec_times_count_fp = cls._spec_times_count_fps[spec_type]
        spec_df = pd.read_csv(spec_fp, sep='\t', encoding="utf-8")
        times_count = spec_df["times"].value_counts()
        sort_times_count = times_count.sort_index()
        cls._write_times_count(sort_times_count, spec_times_count_fp)

    @classmethod
    def process(cls):
        cls._count_spec("cat")
        cls._count_spec("sol")


if __name__ == '__main__':
    SpectatorsCounter.process()
