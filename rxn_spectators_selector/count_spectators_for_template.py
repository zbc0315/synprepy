#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Time    : 2021/9/27 17:09
# @Author  : zhangbc0315@outlook.com
# @File    : count_spectators_for_template.py
# @Software: PyCharm

import pandas as pd

from rmp_database.old_chem_db import OldChemDb
from config import SSPath


class CountSpectatorsForTemplate:

    odb = OldChemDb()

    @classmethod
    def _count_cat_sol_for_tid(cls, df):
        tids = set(df['tid'])
        sids = list(set(df['spectator_id']))
        sorted(sids)
        for tid in tids:
            df_tid = df.query(f'tid=={tid}')
            tot_num = len(df_tid)
            counts = df_tid['spectators'].value_counts()
            tid_sids = [0 if sid not in counts.index else counts[sid] / tot_num for sid in sids]

    @classmethod
    def _get_tid_by_rid(cls, rid):
        rxn_data = cls.odb.search_one_data('public.clean_duplicate_rxn', ['tid_00'], f'rid={rid}')
        tid = rxn_data['tid_00']
        return int(tid) if tid is not None else 0

    @classmethod
    def _get_one_hot(cls, max_idx, idxes):
        res = [0]*(max_idx+1)
        idx_count = idxes.value_counts()
        num = idxes.count()
        for idx in idx_count.index:
            res[idx] = idx_count[idx] / num
        return res

    @classmethod
    def _count_static_data_fro_tids(cls, df):
        tids = sorted(set(df.tid))
        max_sid = df.spectator_id.max()
        res = {'tid': [], 'sids': []}
        for tid in tids:
            sids = df.query(f'tid=={tid}').spectator_id
            tid_one_hot = cls._get_one_hot(max_sid, sids)
            res['tid'].append(tid)
            res['sids'].append(tid_one_hot)
        return res

    @classmethod
    def _count_cat_or_sol(cls, rxn_info_fp: str, s_type: str):
        spec_tid_cor_fps = {'cat': SSPath.CAT_TID_COR_FP,
                            'sol': SSPath.SOL_TID_COR_FP}
        rxn_info_df = pd.read_csv(rxn_info_fp, sep='\t', encoding='utf-8')
        if 'tid' not in rxn_info_df.columns:
            rxn_info_df['tid'] = rxn_info_df['rid'].map(cls._get_tid_by_rid)
            rxn_info_df.to_csv(rxn_info_fp, sep='\t', encoding='utf-8', index=False)
        tid_sids = cls._count_static_data_fro_tids(rxn_info_df)
        tid_sids_df = pd.DataFrame(tid_sids)
        tid_sids_df.to_csv(spec_tid_cor_fps[s_type], sep='\t', encoding='utf-8', index=False)

    @classmethod
    def count(cls):
        cls._count_cat_or_sol(rxn_info_fp=SSPath.CAT_GNN_DATA_INFO_FP, s_type='cat')
        cls._count_cat_or_sol(rxn_info_fp=SSPath.SOL_GNN_DATA_INFO_FP, s_type='sol')


if __name__ == "__main__":
    CountSpectatorsForTemplate.count()
