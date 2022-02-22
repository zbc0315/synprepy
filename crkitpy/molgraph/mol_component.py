# -*- coding: utf-8 -*-
# @Time     : 2020/4/20 16:22
# @Author   : zbc
# @Email    : zhangbc0315@outlook.com
# @FileName : mol_component.py
# @Software : PyCharm


class MolComponent:

    def __init__(self):
        self._molecule = None
        self._copy_hash = self.__hash__()
        self._parent = None

    def set_molecule(self, molecule) -> None:
        self._molecule = molecule

    def get_molecule(self):
        return self._molecule

    def set_copy_hash(self, h: int) -> None:
        self._copy_hash = h

    def get_copy_hash(self) -> int:
        return self._copy_hash

    def copy_equal(self, other) -> bool:
        return self.get_copy_hash() == other.get_copy_hash()

    def set_parent(self, parent) -> None:
        self._parent = parent

    def get_parent(self):
        return self._parent
