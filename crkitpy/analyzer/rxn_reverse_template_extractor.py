# -*- coding: utf-8 -*-
# @Time     : 2020/4/24 9:31
# @Author   : zbc
# @Email    : zhangbc0315@outlook.com
# @FileName : rxn_reverse_template_extractor.py
# @Software : PyCharm

from typing import List, Iterator, Tuple

from crkitpy.rxnhandler.rxn_center_detector import RxnCenterDetector
from crkitpy.molgraph.reaction import Reaction
from crkitpy.molgraph.molecule import Molecule
from crkitpy.molgraph.atom import Atom
from crkitpy.molgraph.bond import Bond
from crkitpy.const import ATOM_PROP


class RxnReverseTemplateExtractor:

    @classmethod
    def extract_rxn_templates(cls, reaction: Reaction) -> Iterator[Molecule]:
        RxnCenterDetector.detect_rxn_center(reaction)
        # cls._detect_dropped_atoms(reaction)
        for rxn_center in cls._extract_rxn_centers(reaction):
            cls._detect_used_linked_bonds(reaction, rxn_center)
            cls._detect_used_broken_bonds(reaction, rxn_center)
            yield cls._extract_rxn_template(reaction, rxn_center)

    # ====================
    # -- 抽取所有反应中心 --

    @classmethod
    def _extract_rxn_centers(cls, reaction: Reaction) -> Iterator[Molecule]:
        """ 抽取反应中心集合
        1. 从产物中抽取；
        2. 依次处理产物中的反应区域（所谓反应区域即为相互关联的反应点（所谓反应点即为发生了断键连键的原子）的集合）；
        3. 每个反应区域向外膨胀一层原子；
        4. 每个反应区域再向外膨胀一层原子，并将这层原子标记为IS RGROUP（即忽视这些原子的元素类型等化学属性）；
        5. 将被标记为 IS RGROUP 的原子的符号设置为‘*’；

        :param reaction: 待抽取反应中心的化学反应
        :return: 反应中心的集合
        """
        for product in reaction.get_products():
            for reacted_atoms in cls._gather_neighbored_reacted_atoms(product, reaction):
                reacted_atoms = cls._expand_atoms(reacted_atoms, 1)
                reacted_atoms = cls._expand_atoms(reacted_atoms, 1, (ATOM_PROP.K_IS_RGROUP, True))
                rxn_center = product.copy_sub_mol(reacted_atoms)
                cls._reset_rgroup_symbol(rxn_center.get_atoms())
                for atom in reacted_atoms:
                    if atom.get_property(ATOM_PROP.K_IS_RGROUP):
                        atom.set_property(ATOM_PROP.K_IS_RGROUP, None)
                yield rxn_center

    @staticmethod
    def _reset_rgroup_symbol(atoms: List[Atom]) -> None:
        """ 将atoms中IS RGROUP的atom的symbol设置为'*'

        :param atoms:
        :return:
        """
        for atom in atoms:
            if atom.get_property(ATOM_PROP.K_IS_RGROUP):
                atom.set_symbol('*')

    @staticmethod
    def _expand_atoms(atoms: List[Atom], expand_num: int, prop: Tuple[str, object] = None) -> List[Atom]:
        """ 处理atoms，将其中每一个atom的临近原子加入到atoms中，循环此操作的次数为expand_num

        :param atoms:
        :param expand_num:
        :param prop: 默认为None, 类型为Tuple[Key, Value]，所有新增加到atoms的原子都会被设置属性，属性名为Key，属性值为Value
        :return:
        """
        for i in range(expand_num):
            new_atoms = atoms.copy()
            for atom in atoms:
                for neighbor_atom in atom.get_neighbors():
                    if neighbor_atom not in atoms:
                        if prop is not None:
                            neighbor_atom.set_property(prop[0], prop[1])
                        new_atoms.append(neighbor_atom)
            atoms = new_atoms
        return atoms

    @classmethod
    def _gather_neighbored_reacted_atoms(cls, molecule: Molecule, reaction: Reaction) -> List[List[Atom]]:
        """ 将化合物中所有IS_REACTED的原子收集出来，相邻的原子聚集为一组，每组即为一个反应中心
        例如，分子中有原子1~10，其中4，7，8，1，3，10都是IS_REACTED，4与7相邻，7与8相邻，1与3相邻，3与10相邻
             则返回结果为[[Atom(4), Atom(7), Atom(8),
                        [Atom(1), Atom(3), Atom(10)]]

        :param molecule: 原子
        :return: 反应中心涉及的原子集合的集合
        """
        reacted_atoms = []
        for atom in molecule.get_atoms():
            if atom.get_property(ATOM_PROP.K_IS_REACTED):
                reacted_atoms.append(atom)
        gather_reacted_atoms = [[atom] for atom in reacted_atoms]
        find_neighbored = True
        while find_neighbored:
            find_neighbored = False
            new_gather_reacted_atoms = []
            gathered_atoms = []
            for i, i_atoms in enumerate(gather_reacted_atoms):
                if i_atoms[0] not in gathered_atoms:
                    new_gather_reacted_atoms.append(i_atoms)
                    index = len(new_gather_reacted_atoms) - 1
                    gathered_atoms.extend(i_atoms)
                else:
                    index = cls._get_index(new_gather_reacted_atoms, i_atoms)
                for j, j_atoms in enumerate(gather_reacted_atoms):
                    if j <= i:
                        continue
                    if cls._is_link(i_atoms, j_atoms, reaction):
                        find_neighbored = True
                        new_gather_reacted_atoms[index].extend(j_atoms)
                        new_gather_reacted_atoms[index] = list(set(new_gather_reacted_atoms[index]))
                        gathered_atoms.extend(j_atoms)
            gather_reacted_atoms = new_gather_reacted_atoms
        return gather_reacted_atoms

    @staticmethod
    def _get_index(gather_reacted_atoms: List[List[Atom]], atoms: List[Atom]) -> int:
        for i, gather_reacted_atom in enumerate(gather_reacted_atoms):
            if atoms[0] in gather_reacted_atom:
                return i

    @classmethod
    def _is_link(cls, i_atoms: List[Atom], j_atoms: List[Atom], reaction: Reaction) -> bool:
        for i_atom in i_atoms:
            for j_atom in j_atoms:
                if i_atom.is_neighbored(j_atom) or i_atom == j_atom:
                    return True
                elif cls._is_linked_in_reactants(i_atom, j_atom, reaction):
                    return True
        return False

    @classmethod
    def _is_linked_in_reactants(cls, i_atom: Atom, j_atom: Atom, reaction: Reaction) -> bool:
        r_i_atom = cls._get_mapping_atom_from_reactants(reaction, i_atom.map_num)
        r_j_atom = cls._get_mapping_atom_from_reactants(reaction, j_atom.map_num)
        if r_i_atom.is_neighbored(r_j_atom) or r_i_atom == r_j_atom:
            return True
        else:
            return False

    # =================================
    # -- 探测反应过程中被舍弃的反应物原子 --

    @staticmethod
    def _detect_dropped_atoms(reaction: Reaction) -> None:
        pass

    # ===============================
    # -- 探测反应过程中被断开的反应物键 --

    @classmethod
    def _detect_used_linked_bonds(cls, reaction: Reaction, rxn_center: Molecule) -> None:
        for c_atom in rxn_center.get_atoms():
            if not c_atom.get_property(ATOM_PROP.K_IS_REACTED):
                continue
            r_atom = cls._get_mapping_atom_from_reactants(reaction, c_atom.map_num)
            for r_n_atom in r_atom.get_neighbors():
                if r_n_atom.map_num is None or r_n_atom.map_num == 0:
                    new_r_n_atom = Atom.instance_by_symbol(r_n_atom.get_symbol())
                    new_r_n_atom.set_property(ATOM_PROP.K_MAYBE_RGROUP, True)
                    rxn_center.add_used_linked_bond_by_type(r_atom.get_bond(r_n_atom).get_bond_type(), c_atom,
                                                            new_r_n_atom)
                else:
                    c_n_atom = cls._get_mapping_atom_from_atoms(c_atom.get_neighbors(), r_n_atom.map_num)
                    r_bond = r_atom.get_bond(r_n_atom)
                    if c_n_atom is None:
                        c_map_n_atoms = list(rxn_center.get_atoms_by_map_num(r_n_atom.map_num))
                        for c_map_n_atom in c_map_n_atoms:
                            # TODO 只添加一次
                            rxn_center.add_used_linked_bond_by_type(r_bond.get_bond_type(), c_atom, c_map_n_atom)
                    else:
                        c_bond = c_atom.get_bond(c_n_atom)
                        if c_bond.get_bond_type() != r_bond.get_bond_type():
                            rxn_center.add_used_linked_bond_by_type(r_bond.get_bond_type(), c_atom, c_n_atom)
        cls._duplicate_and_check_used_linked_bonds(rxn_center)

    @classmethod
    def _duplicate_and_check_used_linked_bonds(cls, rxn_center: Molecule) -> None:
        bonds = rxn_center.get_used_linked_bonds().copy()
        rxn_center.clear_used_linked_bonds()
        added_bonds = []
        for bond in bonds:
            if bond in added_bonds:
                continue
            if any([atom.map_num is None for atom in bond.get_atoms()]):
                rxn_center.add_used_linked_bond(bond)
            else:
                equal_bonds = cls._get_equal_bonds_by_type_and_atoms(bonds, bond)
                if len(equal_bonds) == 2:
                    added_bonds.extend(equal_bonds)
                    rxn_center.add_used_linked_bond(bond)
                elif len(equal_bonds) > 2:
                    raise ValueError('something unknown wrong')

    @staticmethod
    def _get_equal_bonds_by_type_and_atoms(bonds: List[Bond], bond: Bond) -> List[Bond]:
        equal_bonds = []
        for _bond in bonds:
            if bond.equal_in_type_and_atoms(_bond):
                equal_bonds.append(_bond)
        return equal_bonds

    @staticmethod
    def _get_mapping_atom_from_atoms(atoms: Iterator[Atom], map_num: int) -> Atom:
        for atom in atoms:
            if atom.map_num == map_num:
                return atom

    @staticmethod
    def _get_mapping_atom_from_reactants(reaction: Reaction, map_num: int, accept_none=False) -> Atom:
        for reactant in reaction.get_reactants():
            for atom in reactant.get_atoms_by_map_num(map_num):
                if atom is not None:
                    return atom
        if accept_none:
            return None
        else:
            raise ValueError('Cannot find mapping atom from reactants by map num {}'.format(map_num))

    # =============================
    # -- 探测反应过程中被新建连接的键 --

    @classmethod
    def _detect_used_broken_bonds(cls, reaction: Reaction, rxn_center: Molecule) -> None:
        for c_atom in rxn_center.get_atoms():
            if not c_atom.get_property(ATOM_PROP.K_IS_REACTED):
                continue
            for c_n_atom in c_atom.get_neighbors():
                c_bond = c_atom.get_bond(c_n_atom)
                r_atom = cls._get_mapping_atom_from_reactants(reaction, c_atom.map_num)
                r_n_atom = cls._get_mapping_atom_from_reactants(reaction, c_n_atom.map_num, accept_none=True)
                if r_n_atom is None:
                    raise ValueError('Cannot get mapping atom from reactant by product atom')
                r_bond = r_atom.get_bond(r_n_atom)
                if r_bond is None or r_bond.get_bond_type() != c_bond.get_bond_type():
                    if c_bond not in rxn_center.get_used_broken_bonds():
                        rxn_center.add_used_broken_bond(c_bond)

    # =============================
    # -- 根据一个反应中心抽取反应模板 --

    @staticmethod
    def _extract_rxn_template(reaction: Reaction, rxn_center: Molecule) -> Molecule:
        return rxn_center
