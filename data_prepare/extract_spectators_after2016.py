# -*- coding: utf-8 -*-
# @Time     : 2021/4/16 14:59
# @Author   : zbc
# @Email    : zhangbc0315@outlook.com
# @FileName : extract_spectators_after2016.py
# @Software : PyCharm

from typing import Iterator
import json

from rmp_database.old_chem_db import OldChemDb
from log_utils import Logger, Note
from config import SSPath


class ExtractSpectatorsAfter2016:
    """ 处理2017年及以后的数据，将催化剂、溶剂信息提取出来

    """
    _odb = OldChemDb()
    _spectators = {"sid": [], "inchi": [], "smiles": [], "been_catalyst_times": [], "been_solvent_times": []}
    _new_spectators_num = 0

    @classmethod
    def _get_all_saved_spectators(cls):
        """ 从 public.spectators_mols 中获得所有催化剂、溶剂

        :return:
        """
        n = 0
        for spectators_data in cls._odb.get_data_iter("public.spectators_mols",
                                                      ["sid", "inchi", "smiles",
                                                       "been_catalyst_times", "been_solvent_times"],
                                                      None):
            n += 1
            for key in spectators_data.keys():
                cls._spectators[key].append(spectators_data[key])
        Note.note(SSPath.SS_NOTE_FP, f"Spectators Before 2016: {n}")
        Logger.info(f"Load Spectators: {n}")

    @classmethod
    def _get_rids(cls) -> Iterator[int]:
        """ 从 rxn_for_test.clean_duplicate_rxn 中获得所有的反应的 rid

        :return: Iterator[rid]
        """
        odb = OldChemDb()
        for rxn_data in odb.get_data_iter('rxn_for_test.clean_duplicate_rxn', ['rid'], None):
            yield rxn_data['rid']

    @classmethod
    def _parse_spectators(cls, spectators_str: str) -> [(str, str, str)]:
        """ 解析spectators字符串

        :param spectators_str:
        :return: [(name, smiles, inchi)]
        """
        if spectators_str is None or len(spectators_str) == 0:
            return []
        res = []
        spectators = json.loads(spectators_str)
        for spectator in spectators:
            if "smiles" not in spectator.keys():
                continue
            res.append((spectator["name"], spectator["smiles"], spectator['inchi']))
        return res

    @classmethod
    def _save_spectator(cls, mol: (str, str, str), role: str) -> int:
        """ 如果inchi无法在cls._spectators中查到，则存储 催化剂、溶剂到 public.spectators_mols

        :param mol: (name, smiles, inchi)
        :return: 分子的sid
        """
        been_catalyst_times = 0
        been_solvent_times = 0
        if role == "cat":
            been_catalyst_times += 1
        elif role == "sol":
            been_solvent_times += 1
        else:
            raise ValueError(f"Unexpected role: {role}")
        if mol[2] not in cls._spectators['inchi']:
            Logger.info("new mol")
            cls._new_spectators_num += 1
            sid = len(cls._spectators['sid']) + 1
            cls._odb.insert_one_data("public.spectators_mols",
                                     ["sid", "name", "smiles", "inchi", "been_catalyst_times", "been_solvent_times"],
                                     [sid, mol[0], mol[1], mol[2], been_catalyst_times, been_solvent_times])
            cls._spectators['sid'].append(sid)
            cls._spectators['smiles'].append(mol[1])
            cls._spectators['inchi'].append(mol[2])
            cls._spectators['been_catalyst_times'].append(been_catalyst_times)
            cls._spectators['been_solvent_times'].append(been_solvent_times)
        else:
            Logger.info("mol exists")
            idx = cls._spectators['inchi'].index(mol[2])
            sid = cls._spectators['sid'][idx]
            cls._spectators['been_catalyst_times'][idx] += been_catalyst_times
            cls._spectators['been_solvent_times'][idx] += been_solvent_times
        return sid

    @classmethod
    def _save_spectators(cls, mols: [(str, str, str)], role: str) -> [int]:
        """ 如果inchi无法在cls._spectators中查到，则存储 催化剂、溶剂到 public.spectators_mols

        :param mols:
        :param role:
        :return: 所有传入分子的sid列表
        """
        res = []
        for mol in mols:
            res.append(cls._save_spectator(mol, role))
        return res

    @classmethod
    def _get_spectators_by_rid(cls, rid: int) -> ([(str, str, str)], [(str, str, str)]):
        """ 从 rxn_for_test.origin_rxn 中获得 rid 对应的反应的催化剂、溶剂信息

        :param rid:
        :return: [(cat_name, cat_smiles, cat_inchi)], [(sol_name, sol_smiles, sol_inchi)]
        """
        rxn_data = cls._odb.search_one_data("rxn_for_test.origin_rxn", ["catalysts", "solvents"], f"rid={rid}")
        return cls._parse_spectators(rxn_data["catalysts"]), cls._parse_spectators(rxn_data["solvents"])

    @classmethod
    def _get_mols_code(cls, ids: [int]) -> str:
        """

        :param ids: [1,2,3]
        :return: "1.2.3"
        """
        if len(ids) == 0:
            return ''
        return '.'.join([str(i) for i in ids])

    @classmethod
    def process(cls):
        cls._get_all_saved_spectators()
        for rid in cls._get_rids():
            cats, sols = cls._get_spectators_by_rid(rid)
            cats_sids = cls._save_spectators(cats, 'cat')
            sols_sids = cls._save_spectators(sols, 'sol')

            cats_sids = list(set(cats_sids))
            sols_sids = list(set(sols_sids))
            if len(cats_sids) > 0:
                print(cats_sids)
            cls._odb.update("rxn_for_test.clean_duplicate_rxn",
                            ["cats", "sols"],
                            [cls._get_mols_code(cats_sids), cls._get_mols_code(sols_sids)],
                            f"rid={rid}")
        Note.note(SSPath.SS_NOTE_FP, f"Spectators After 2016 More: {cls._new_spectators_num}")
        for i in range(len(cls._spectators['sid'])):
            sid = cls._spectators['sid'][i]
            cat_times = cls._spectators['been_catalyst_times'][i]
            sol_times = cls._spectators['been_solvent_times'][i]
            cls._odb.update("public.spectators_mols",
                            ['been_catalyst_times', 'been_solvent_times'],
                            [cat_times, sol_times],
                            f'sid={sid}')


if __name__ == '__main__':
    ExtractSpectatorsAfter2016.process()
