# -*- coding: utf-8 -*-
# @Time     : 2020/5/28 16:57
# @Author   : zbc
# @Email    : zhangbc0315@outlook.com
# @FileName : rsmates_parser.py
# @Software : PyCharm

from typing import Tuple

from crkitpy.smilesparse.smiles_parser import SmilesParser
from crkitpy.molgraph.molecule import Molecule, Bond, Atom
from crkitpy.const import BOND_TYPE


class RSmatesParser:
    """ Smates是本程序定义的，类似于smiles，是template的唯一字符串表示

    """

    @classmethod
    def smates_to_template(cls, smates: str):
        smiles, \
            invisible_info, \
            used_linked_bond_info, \
            used_broken_bond_info, \
            changed_charge_info, \
            changed_valence_info = smates.split(' | ')
        template = SmilesParser.parse_to_mol_by_lex(smiles)
        cls._parse_invisible_info(template, invisible_info)
        cls._parse_used_linked_bond_info(template, used_linked_bond_info)
        cls._parse_used_broken_bond_info(template, used_broken_bond_info)
        cls._parse_changed_charge_info(template, changed_charge_info)
        cls._parse_changed_valence_info(template, changed_valence_info)
        return template

    # ================================
    # -- parse changed valence info --

    @classmethod
    def _parse_changed_valence_info(cls, template: Molecule, changed_valence_info: str):
        changed_valence_info = changed_valence_info.lstrip('cv:')
        if len(changed_valence_info) == 0:
            return
        for changed_valence in changed_valence_info.split(','):
            map_num, old_valence, new_valence = cls._parse_prop(changed_valence)
            for atom in template.get_atoms_by_map_num(map_num):
                template.add_changed_valence(atom, int(old_valence), int(new_valence))
                break

    # ===============================
    # -- parse changed charge info --

    @classmethod
    def _parse_changed_charge_info(cls, template: Molecule, changed_charge_info: str):
        changed_charge_info = changed_charge_info.lstrip('cc:')
        if len(changed_charge_info) == 0:
            return
        for changed_charge in changed_charge_info.split(','):
            map_num, old_charge, new_charge = cls._parse_prop(changed_charge)
            for atom in template.get_atoms_by_map_num(map_num):
                template.add_changed_charge(atom, int(old_charge), int(new_charge))
                break

    @staticmethod
    def _parse_prop(prop_info: str) -> Tuple[int, int, int]:
        if '--' in prop_info:
            prop_info = prop_info.replace('-', '=')
            prop_info = prop_info.replace('==', '=-')
            i, j, k = prop_info.split('=')
        else:
            i, j, k = prop_info.split('-')
        return int(i), int(j), int(k)

    # ========================================
    # -- parse used linked/broken bond info --

    @classmethod
    def _parse_used_linked_bond_info(cls, template: Molecule, used_linked_bond_info: str):
        used_linked_bond_info = used_linked_bond_info.lstrip('ul:')
        if len(used_linked_bond_info) == 0:
            return
        for used_linked_bond in used_linked_bond_info.split(','):
            begin_atom, bond_type, end_atom = cls._parse_bond(used_linked_bond, template)
            template.add_used_linked_bond_by_type(bond_type, begin_atom, end_atom)

    @classmethod
    def _parse_used_broken_bond_info(cls, template: Molecule, used_broken_bond_info: str):
        used_broken_bond_info = used_broken_bond_info.lstrip('ub:')
        if len(used_broken_bond_info) == 0:
            return
        for used_broken_bond in used_broken_bond_info.split(','):
            begin_atom, bond_type, end_atom = cls._parse_bond(used_broken_bond, template)
            bond = template.get_bond(begin_atom, end_atom)
            template.add_used_broken_bond(bond)

    @staticmethod
    def _parse_bond(bond_info: str, template: Molecule) -> (Atom, str, Atom):
        begin_map_num, bond_weight, end_map_num = bond_info.split('-')
        begin_atom = list(template.get_atoms_by_map_num(int(begin_map_num)))[0]
        end_atom = list(template.get_atoms_by_map_num(int(end_map_num)))[0]
        bond_type = BOND_TYPE.BOND_WEIGHT_STR_TO_SYMBOLS[bond_weight]
        return begin_atom, bond_type, end_atom

    # ==========================
    # -- parse invisible info --

    @classmethod
    def _parse_invisible_info(cls, template: Molecule, invisible_info: str):
        invisible_info = invisible_info.lstrip('iv:')
        if len(invisible_info) == 0:
            return
        map_nums = [int(map_num) for map_num in invisible_info.split(',')]
        for atom in template.get_atoms_by_map_nums(map_nums):
            atom.is_visible = False
