#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Time    : 2022/1/5 13:12
# @Author  : zhangbc0315@outlook.com
# @File    : atom_handler.py
# @Software: PyCharm

from rdkit.Chem.rdchem import Atom


class AtomHandler:

    TYPE_DROPPED = "dropped"
    TYPE_REACTED = "reacted"
    TYPE_NBR = "nbr"

    @classmethod
    def set_atom_is_nbr(cls, atom: Atom):
        atom.SetProp("type", cls.TYPE_NBR)

    @classmethod
    def is_atom_nbr(cls, atom: Atom) -> bool:
        return "type" in atom.GetPropNames() and atom.GetProp("type") == cls.TYPE_NBR

    @classmethod
    def set_atom_is_reacted(cls, atom: Atom):
        atom.SetProp("type", cls.TYPE_REACTED)

    @classmethod
    def is_atom_reacted(cls, atom: Atom) -> bool:
        return "type" in atom.GetPropNames() and atom.GetProp("type") == cls.TYPE_REACTED

    @classmethod
    def set_atom_is_dropped(cls, atom: Atom):
        atom.SetProp("type", cls.TYPE_DROPPED)

    @classmethod
    def is_atom_dropped(cls, atom: Atom) -> bool:
        return "type" in atom.GetPropNames() and atom.GetProp("type") == cls.TYPE_DROPPED


if __name__ == "__main__":
    pass
