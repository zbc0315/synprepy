# -*- coding: utf-8 -*-
# @Time     : 2020/7/1 14:52
# @Author   : zbc
# @Email    : zhangbc0315@outlook.com
# @FileName : electron_state.py
# @Software : PyCharm


class ElectronState:

    _l_str = ['s', 'p', 'd', 'f', 'g', 'h', 'i', 'j', 'k']

    def __init__(self, n: int = None, l: int = None, m: int = None, ms: float = None, num: int = None):
        self._n   = n
        self._l   = l
        self._m   = m
        self._ms  = ms
        self._num = num

    def __str__(self):
        state_str = None
        if self.get_l() is not None:
            state_str = f'{self.get_n()}{self._l_str[self.get_l()]}'
        elif self.get_n() is not None:
            state_str = str(self.get_n())
        if self.get_num() is not None:
            return f'{state_str}({self.get_num()})'
        else:
            return state_str

    def set_n(self, n: int):
        self._n = n

    def get_n(self) -> int:
        return self._n

    def set_l(self, l: int):
        self._l = l

    def get_l(self) -> int:
        return self._l

    def set_m(self, m: int):
        self._m = m

    def get_m(self) -> int:
        return self._m

    def set_ms(self, ms: float):
        self._ms = ms

    def get_ms(self) -> float:
        return self._ms

    def set_num(self, num: int):
        self._num = num

    def get_num(self) -> int:
        return self._num

    def get_max_num(self) -> int:
        if self.get_ms() is not None:
            return 1
        elif self.get_m() is not None:
            return 2
        elif self.get_l() is not None:
            return (self.get_l() + 1) * 2
        elif self.get_n() is not None:
            return (self.get_n() ** 2) * 2
