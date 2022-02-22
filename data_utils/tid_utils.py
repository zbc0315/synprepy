# -*- coding: utf-8 -*-
# @Time     : 2020/11/9 15:54
# @Author   : zbc
# @Email    : zhangbc0315@outlook.com
# @FileName : tid_utils.py
# @Software : PyCharm

from typing import Union, List, Dict

from config import BTPath
from config import CTPath


class TidUtils:
    _tids_list: Union[None, List[int]] = None
    _tids_dict: Dict[int, int] = {}

    @classmethod
    def load_comp_tid_to_idx(cls, fp: str = None, res_type: str = 'dict'):
        if fp is None:
            return cls._load_tids(CTPath.COMP_TIDS_FP, sep='\n', res_type=res_type)
        else:
            return cls._load_tids(fp, sep='\n', res_type=res_type)

    @classmethod
    def save_comp_tids(cls, tids):
        cls._save_tids(tids, CTPath.COMP_TIDS_FP)

    # ==================================

    @classmethod
    def load_basic_tids(cls, fp: str = None, sep: str = '.'):
        if fp is None:
            return cls._load_tids(BTPath.TIDS_FP, sep=sep)
        else:
            return cls._load_tids(fp, sep=sep)

    @classmethod
    def load_tids(cls):
        return cls._load_tids(BTPath.TIDS_FP)

    # ===================================

    @classmethod
    def _save_tids(cls, tids, fp):
        tids_str = '\n'.join([str(tid) for tid in tids])
        with open(fp, 'w', encoding='utf-8')as f:
            f.write(tids_str)

    @classmethod
    def _load_tids(cls, fp: str, sep: str = '.', res_type: str = 'list'):
        cls._tids_list = []
        with open(fp, 'r', encoding='utf-8')as f:
            for n, tid in enumerate(f.read().split(sep)):
                cls._tids_list.append(int(tid))
                cls._tids_dict[int(tid)] = n
        if res_type == 'list':
            return cls._tids_list
        else:
            return cls._tids_dict

    @classmethod
    def encode_tids(cls, tids: List[int]) -> List[float]:
        if cls._tids_list is None:
            cls.load_tids()
        prop = 1 / len(tids)
        res = [0.0] * len(cls._tids_list)
        for tid in tids:
            idx = cls._tids_dict[tid]
            res[idx] = prop
        return res


if __name__ == '__main__':
    pass
