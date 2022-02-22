# -*- coding: utf-8 -*-
# @Time     : 2020/5/28 16:58
# @Author   : zbc
# @Email    : zhangbc0315@outlook.com
# @FileName : rsmates_writer.py
# @Software : PyCharm

from typing import List, Dict, Iterator

from crkitpy.smilesparse.smiles_writer import SmilesWriter
from crkitpy.molgraph.molecule import Molecule, Atom, Bond
from crkitpy.const import ATOM_PROP, BOND_TYPE
from crkitpy.rxnhandler.map_num_reset import MapNumReset


class RSmatesWriter:

    """ Smates是本程序定义的，类似于smiles，是template的唯一字符串表示
        TODO 未处理分子中的氢
    """

    @classmethod
    def template_to_smates(cls, template: Molecule) -> str:
        MapNumReset.reset_molecule_map_num_to_id(template)
        components = template.get_detached_sub_mols()
        components = cls._sort_components(components)
        smates_sb = [cls._write_smiles(components),
                     cls._write_invisible_info(components),
                     cls._write_used_linked_bond_info(template),
                     cls._write_used_broken_bond_info(template),
                     cls._write_changed_charge_info(template),
                     cls._write_changed_valence_info(template)]
        return ' | '.join(smates_sb)

    # ===========================
    # -- write changed valence --

    @classmethod
    def _write_changed_valence_info(cls, template: Molecule) -> str:
        changed_valence_info_sb = []
        for atom, (old_valence, new_valence) in template.get_changed_valence().items():
            changed_valence_info_sb.append(f'{atom.map_num}-{old_valence}-{new_valence}')
        return 'cv:' + ','.join(changed_valence_info_sb)

    # ==========================
    # -- write changed charge --

    @classmethod
    def _write_changed_charge_info(cls, template: Molecule) -> str:
        changed_charge_info_sb = []
        for atom, (old_charge, new_charge) in template.get_changed_charge().items():
            changed_charge_info_sb.append(f'{atom.map_num}-{old_charge}-{new_charge}')
        return 'cc:' + ','.join(changed_charge_info_sb)

    # ===========================
    # -- write used linked/broken bond info --

    @classmethod
    def _write_used_linked_bond_info(cls, template: Molecule) -> str:
        used_linked_bond_info_sb = []
        for bond in template.get_used_linked_bonds():
            used_linked_bond_info_sb.append(cls._bond_info(bond))
        return 'ul:' + ','.join(used_linked_bond_info_sb)

    @classmethod
    def _write_used_broken_bond_info(cls, tempalte: Molecule) -> str:
        used_broken_bond_info_sb = []
        for bond in tempalte.get_used_broken_bonds():
            used_broken_bond_info_sb.append(cls._bond_info(bond))
        return 'ub:' + ','.join(used_broken_bond_info_sb)

    @staticmethod
    def _bond_info(bond: Bond):
        bond_weight = BOND_TYPE.BOND_WEIGHTS[bond.get_bond_type()]
        map_nums = [a.map_num for a in bond.get_atoms()]
        map_nums.sort()
        return f'{map_nums[0]}-{bond_weight}-{map_nums[1]}'

    # ========================
    # -- write invisible info --

    @classmethod
    def _write_invisible_info(cls, components: List[Molecule]) -> str:
        visible_map_nums = []
        for n, component in enumerate(components):
            if not cls._is_component_visible(component):
                visible_map_nums.extend([atom.map_num for atom in component.get_atoms()])
        visible_map_nums.sort()
        visible_info_sb = [str(v) for v in visible_map_nums]
        return 'iv:'+','.join(visible_info_sb)

    # ==================
    # -- write smiles --

    @classmethod
    def _write_smiles(cls, components: List[Molecule]) -> str:
        smileses = []
        for component in components:
            smileses.append(SmilesWriter.mol_to_smiles(component, need_reset_map_num=False, need_normal_smiles=False))
        return '.'.join(smileses)

    # ==========================
    # -- is component visible --

    @classmethod
    def _is_component_visible(cls, component: Molecule) -> bool:
        if all([atom.is_visible for atom in component.get_atoms()]):
            return True
        elif all([not atom.is_visible for atom in component.get_atoms()]):
            return False
        else:
            raise ValueError('there are visible and invisible atoms exist in one component')

    # =====================
    # -- sort components --

    @classmethod
    def _sort_components(cls, components: Iterator[Molecule]) -> List[Molecule]:
        return sorted(components, key=lambda x: cls._get_min_id(x))

    @classmethod
    def _get_min_id(cls, molecule: Molecule) -> int:
        min_id = None
        for atom in molecule.get_atoms():
            if min_id is None or min_id > atom.get_id():
                min_id = atom.get_id()
        return min_id

