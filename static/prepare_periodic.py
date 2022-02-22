# -*- coding: utf-8 -*-
# @Time     : 2021/3/1 15:49
# @Author   : zbc
# @Email    : zhangbc0315@outlook.com
# @FileName : prepare_periodic.py
# @Software : PyCharm

import json
import pandas as pd


class PreparePeriodic:

    ls = ['s', 'p', 'd', 'f']

    @classmethod
    def parse_nln(cls, nln: str):
        n = int(nln[0])
        l = cls.ls.index(nln[1])
        num = int(nln[2:])
        return [n, l, num]

    @classmethod
    def parse_simple_nln(cls, df: pd.DataFrame, sim_nln: str):
        q = df[df.symbol == sim_nln]
        q = q.reset_index()
        return json.loads(q.loc[0]['configuration'])

    @classmethod
    def parse_configuration(cls, df: pd.DataFrame, conf: str):
        conf = conf.split(' ')
        configuration = []
        for con in conf:
            if con.startswith('['):
                configuration = cls.parse_simple_nln(df, con[1:-1])
            else:
                configuration.append(cls.parse_nln(con))
        return json.dumps(configuration)

    @classmethod
    def get_configuration(cls, df: pd.DataFrame):
        df['configuration'] = [None]*len(df.index)
        for i in df.index:
            df.loc[i, 'configuration'] = cls.parse_configuration(df, df.loc[i]['simple_configuration'])
        return df

    @classmethod
    def process(cls):
        df = pd.read_csv('periodic.tsv', sep='\t', encoding='utf-8')
        df = cls.get_configuration(df)
        df.to_csv('periodic.tsv', sep='\t', encoding='utf-8', index=False)


if __name__ == '__main__':
    PreparePeriodic.process()
