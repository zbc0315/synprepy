# -*- coding: utf-8 -*-
# @Time     : 2020/6/30 11:18
# @Author   : zbc
# @Email    : zhangbc0315@outlook.com
# @FileName : electron_state_utils.py
# @Software : PyCharm

from typing import Iterator, Tuple, List

from crkitpy.electronstate.electron_state import ElectronState


class ElectronStateUtils:

    l_str = ['s', 'p', 'd', 'f', 'g', 'h', 'i', 'l', 'm']

    # ===============================
    # -- get state by electron num --

    @classmethod
    def get_state_by_electron_num(cls, electron_num: int):
        states = cls.get_sorted_nl_states(9)
        # states = [ElectronState(n=nl[0], l=nl[1]) for nl in states]
        nl_nums = []
        for nl in states:
            if electron_num <= 0:
                break
            max_num = cls.get_m_num_by_l(nl[1])
            if electron_num >= max_num:
                num = max_num
            else:
                num = electron_num
            electron_num -= num
            nl_nums.append((nl[0], nl[1], num))
        nl_nums.sort(key=lambda x: x[0])
        states = [ElectronState(n=nl_num[0], l=nl_num[1], num=nl_num[2]) for nl_num in nl_nums]
        return states

    # ===============
    # -- state num --

    @classmethod
    def get_l_num(cls, n: int) -> int:
        return n

    @classmethod
    def get_m_num_by_l(cls, l: int) -> int:
        return (l*2+1)*2

    # ================
    # -- state iter --

    @classmethod
    def get_n_iter(cls) -> Iterator[int]:
        n = 0
        while True:
            n += 1
            yield n

    @classmethod
    def get_l_iter(cls, n: int) -> Iterator[int]:
        for l in range(n):
            yield l

    @classmethod
    def get_m_iter(cls, l) -> Iterator[int]:
        for m in range(-l, l + 1):
            yield m

    @classmethod
    def get_ms_iter(cls) -> Iterator[int]:
        for ms in [-0.5, 0.5]:
            yield ms

    @classmethod
    def get_nlm_state_iter(cls) -> Iterator[Tuple[int, int, int]]:
        for n in cls.get_n_iter():
            for l in cls.get_l_iter(n):
                for m in cls.get_m_iter(l):
                    yield n, l, m

    @classmethod
    def get_nl_state_iter(cls) -> Iterator[Tuple[int, int]]:
        for n in cls.get_n_iter():
            for l in cls.get_l_iter(n):
                yield n, l

    @classmethod
    def get_sorted_nl_states(cls, max_n: int) -> List[Tuple[int, int]]:
        states = []
        for n, l in cls.get_nl_state_iter():
            if n > max_n:
                break
            states.append((n, l))
        states.sort(key=lambda x: (x[0]+x[1], x[0]))
        return states

    @classmethod
    def states_nl_to_str(cls, nl: Tuple[int, int]) -> str:
        return f'{nl[0]}{cls.l_str[nl[1]]}'


if __name__ == '__main__':
    for n in range(173):
        ss = ElectronStateUtils.get_state_by_electron_num(n)
        s_str = [str(s) for s in ss]
        space_num = 4-len(str(n))
        print(f"{n}{' '*space_num}- {' '.join(s_str)}")
