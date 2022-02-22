# -*- coding: utf-8 -*-
# @Time     : 2021/3/3 13:58
# @Author   : zbc
# @Email    : zhangbc0315@outlook.com
# @FileName : spectators_prepare.py
# @Software : PyCharm

import pandas as pd

from rmp_database.old_chem_db import OldChemDb

from data_prepare import origin_data_path as ORIGIN_DATA_PATH


class SpectatorsPrepare:

    @classmethod
    def _parse_known_spectators_file(cls, fp: str):
        inchis = []
        with open(fp, 'r', encoding='utf-8')as f:
            for line in f.readlines():
                line = line.strip()
                if line.startswith('#'):
                    continue
                inchis.append(line.split('\t')[0])
        return inchis

    @classmethod
    def _load_known_spectators_inchis(cls):
        cats = cls._parse_known_spectators_file(ORIGIN_DATA_PATH.KNOWN_CAT_INCHIS)
        sols = cls._parse_known_spectators_file(ORIGIN_DATA_PATH.KNOWN_SOL_INCHIS)
        return cats, sols

    @classmethod
    def _is_need_change(cls, inchi: str,
                        been_cat_times: int, been_sol_times: int,
                        cat_inchis: [str], sol_inchis: [str]):
        if been_cat_times > been_sol_times and been_sol_times <= 2 and inchi in cat_inchis:
            return 'catalyst'
        elif been_sol_times > been_cat_times and been_cat_times <= 2 and inchi in sol_inchis:
            return 'solvent'
        else:
            return None

    @classmethod
    def clean_spectators(cls):
        cat_inchis, sol_inchis = cls._load_known_spectators_inchis()

        db = OldChemDb()
        for spec_data in db.get_data_iter("public.spectators_mols",
                                          ['inchi', 'been_catalyst_times', 'been_solvent_times'], None):
            been_cat_times = spec_data['been_catalyst_times']
            been_sol_times = spec_data['been_solvent_times']
            if been_cat_times != 0 and been_sol_times != 0:
                spec_type = cls._is_need_change(spec_data['inchi'],
                                                been_cat_times, been_sol_times,
                                                cat_inchis, sol_inchis)

                if spec_type is None:
                    inchi = spec_data['inchi']
                    print(f"{been_cat_times} - {been_sol_times} - {inchi in cat_inchis} - {inchi in sol_inchis}")
                else:
                    print('*'*30)

    @classmethod
    def get_all_rxn_with_spectators(cls):
        all_rxn_df = pd.read_csv(ORIGIN_DATA_PATH.ALL_RXN_FP, sep='\t', encoding='utf-8')
        print(all_rxn_df.size)
        all_rxn_with_spec = all_rxn_df[pd.notna(all_rxn_df.catalysts_sids) | pd.notna(all_rxn_df.solvents_sids)]
        print(all_rxn_with_spec.size)
        all_rxn_with_spec.to_csv(ORIGIN_DATA_PATH.ALL_RXN_WITH_SPEC_FP, sep='\t', encoding='utf-8', index=False)


if __name__ == '__main__':
    SpectatorsPrepare.clean_spectators()
    SpectatorsPrepare.get_all_rxn_with_spectators()
