# -*- coding: utf-8 -*-
# @Time     : 2020/4/20 16:22
# @Author   : zbc
# @Email    : zhangbc0315@outlook.com
# @FileName : atom.py
# @Software : PyCharm

from typing import Iterator, List, Any

from enum import Enum

from crkitpy.molgraph.mol_component import MolComponent
from crkitpy.molgraph.element import Element
from crkitpy.const import BOND_TYPE
from crkitpy.errors.atom_error import *


class Hybrid(Enum):
    UNSPECIFIC = 0
    S = 1
    SP = 2
    SP2 = 3
    SP3 = 4
    SP3D = 5
    SP3D2 = 6


class Atom(MolComponent):
    """ 原子

    """

    def __init__(self):
        super().__init__()
        self._id = None
        self._atomic_num = None
        self._hybrid = Hybrid.UNSPECIFIC
        self._alter_hybrids: List[Hybrid] = []
        self._bonds = []
        self.charge: int = 0
        self.explicit_hs_num = 0
        self.implicit_hs_num = 0
        self._explicit_valence = 0
        self._implicit_valence = 0  # 此处应根据原子所连H原子数决定，目前初始化为0，实际应初始化为None（在解析smiles时彻底替代换下rdkit之后）
        self.is_aromatic = False
        self._symbol = None
        self.map_num = None
        self.degree = None
        self._aromatic_cycles = []
        self.properties = {}
        self._brock_bonds = []
        self.is_visible = True
        self._copy_hash = self.__hash__()

    def __str__(self):
        atom_str = self.get_symbol()
        if self.charge is not None:
            if self.charge > 0:
                atom_str = f'{atom_str}+{self.charge}'
            elif self.charge < 0:
                atom_str = f'{atom_str}{self.charge}'
        if self.map_num is not None:
            atom_str = f'{atom_str}:{self.map_num}'
        return atom_str

    # ========
    # -- id --

    def set_id(self, i: int) -> None:
        self._id = i

    def get_id(self) -> int:
        return self._id

    # ============
    # -- hybrid --

    def set_hybrid(self, hybrid: Hybrid) -> None:
        self._hybrid = hybrid

    def get_hybrid(self) -> Hybrid:
        return self._hybrid

    # =========================
    # -- alternative hybrids --

    def add_alter_hybrid(self, hybrid: Hybrid) -> None:
        if hybrid not in self._alter_hybrids:
            self._alter_hybrids.append(hybrid)

    def get_alter_hybrids(self) -> List[Hybrid]:
        return self._alter_hybrids

    # ===========================
    # -- Symbol and Atomic Num --

    def set_symbol(self, symbol: str) -> None:
        """ 设置元素符号，会同时设置元素序数

        :param symbol: 元素符号，首字母大写; '*'代表通用元素
        """
        self._symbol = symbol
        if self.is_r_group():
            self.set_atomic_num(0)
            self.charge = 0
            self.implicit_hs_num = 0
            self.explicit_hs_num = 0
            self._explicit_valence = 0
            self._implicit_valence = 0
        else:
            self.set_atomic_num(Element.get_atomic_num_by_symbol(symbol))

    def get_symbol(self) -> str:
        return self._symbol

    def is_r_group(self) -> bool:
        return self.get_symbol() in ['*', 'ARO']

    def set_atomic_num(self, atomic_num: int) -> None:
        self._atomic_num = atomic_num

    def get_atomic_num(self) -> int:
        return self._atomic_num

    # ====================
    # -- aromatic cycle --

    def _add_aromatic_cycle(self, cycle) -> None:
        if cycle not in self._aromatic_cycles:
            self._aromatic_cycles.append(cycle)

    def get_aromatic_cycles(self) -> List:
        """ 获得芳环的序列

        :return: List[AtomSet]，芳香环的序列
        """
        return self._aromatic_cycles

    # ======================
    # -- explicit valence --

    def set_explicit_valence(self, ev: int) -> None:
        """ 设置显式价态

        :param ev: 显式价态
        """
        self._explicit_valence = ev

    def get_explicit_valence(self) -> int:
        """ 获得显式价态

        :return: 显式价态
        """
        return self._explicit_valence

    # ======================
    # -- implicit valence --

    def set_implicit_valence(self, iv: int) -> None:
        """ 设置隐式价态

        :param iv: 隐式价态
        """
        self._implicit_valence = iv

    def get_implicit_valence(self) -> int:
        """ 获得隐式价态

        :return: 隐式价态
        """
        return self._implicit_valence

    def get_valence(self) -> int:
        """ 获得价态
        价态 = 隐式价态 + 显式价态

        :return: 价态
        """
        return self.get_explicit_valence() + self.get_implicit_valence()

    # ===============
    # -- Hs number --

    def get_hs_num(self):
        """ 获得氢数
        氢数 = 隐式氢数 + 显式氢数

        :return: 氢数
        """
        if self.implicit_hs_num is not None and self.explicit_hs_num is not None:
            return self.implicit_hs_num + self.explicit_hs_num
        else:
            return 0

    def correct_implicit_hs_num(self):
        """ 纠正隐式氢数
        在分子进行断/连键后，对相应原子添加/删除氢原子，达到维持该原子价态不变的目的
        隐式氢数 = 价态 - 显式氢数 - 总非氢键级
        如果隐式氢数<0，则调整显式氢数，使得隐式氢数为0; 如果显式氢数因此而<0，则将显式氢数设置为0

        调整价态，将隐式价态设置为氢数

        """
        bonds_weight = self.get_bonds_weight_without_h()
        self.implicit_hs_num = self.get_valence() - self.explicit_hs_num - bonds_weight
        if self.implicit_hs_num < 0:
            self.explicit_hs_num += self.implicit_hs_num
            self.implicit_hs_num = 0
        if self.explicit_hs_num < 0:
            self.explicit_hs_num = 0  # 此处直接忽略了错误，可能会导致其他错误
            # raise AtomValenceError('implicit H num can not be negative')
        valence = self.get_implicit_valence() + self.get_explicit_valence()
        self.set_implicit_valence(self.explicit_hs_num + self.implicit_hs_num)
        self.set_explicit_valence(valence - self.get_implicit_valence())

    # ===========
    # -- Envs --

    def get_map_num_env(self):
        """ 获得原子的环境map_num

        :return: 本原子临近原子的map_num有序序列
        """
        neighbor_map_nums = [neighbor.map_num for neighbor in self.get_neighbors()]
        neighbor_map_nums.sort()
        return neighbor_map_nums

    def get_symbol_env(self):
        """ 得到临近的原子的symbol序列
        # TODO 并未考虑临近的键的类型，在用作子图同构判断时或许有用
        :return:
        """
        neighbor_symbols = [neighbor.get_symbol() for neighbor in self.get_neighbors()]
        neighbor_symbols.sort()
        return neighbor_symbols

    # ==========
    # -- Bond --

    def add_bond(self, bond) -> None:
        """ 原子添加化学键，并且化学键添加原子

        :param bond:
        :return:
        """
        self._bonds.append(bond)
        bond.add_atom(self)

    def get_bonds(self) -> List:
        return self._bonds

    def get_bonds_weight_without_h(self) -> int:
        """ 获得除了氢之外的总键级
        芳香键的键级为1.5
        如果是芳香C原子，键级-1
        有奇数个芳香键时键级-0.5，有偶数个芳香键时键级不变
        如果有两个芳香键和一个双键或三键，键级-1

        :return: 除了氢之外的总键级
        """
        bonds_weight = 0
        aromatic_num = 0
        un_single_num = 0
        for bond in self.get_bonds():
            if any([atom.get_symbol() == 'H' for atom in bond.get_atoms()]):
                continue
            else:
                bond_weight = BOND_TYPE.BOND_WEIGHTS[bond.get_bond_type()]
                bonds_weight += bond_weight
                if bond_weight == 1.5:
                    aromatic_num += 1
                elif bond_weight != 1:
                    un_single_num += 1
        if self.get_symbol() != 'C' and self.is_aromatic:
            bonds_weight -= 1
        bonds_weight -= (aromatic_num%2)*0.5
        if aromatic_num == 2 and un_single_num == 1:
            bonds_weight -= 1
        if int(bonds_weight) != bonds_weight:
            raise AtomValenceError('bonds weight can not be decimal')
        return int(bonds_weight)

    def get_bond(self, other_atom):
        """ 获得本原子与其他原子之间的键

        :param other_atom: 其他原子
        :return:
        """
        for bond in self.get_bonds():
            if other_atom.is_linked(bond):
                return bond

    def is_linked(self, bond) -> bool:
        """ 判断本原子是否与键相连

        :param bond: 键
        :return: 是/否相连
        """
        for _bond in self.get_bonds():
            if _bond == bond:
                return True
        return False

    def link_bond(self, bond, other_atom):
        """ 通过molecule连接bond，bond会被加入molecule

        :param bond: 键
        :param other_atom: 本原子通过键连接的另一个原子，other_atom应该存在于molecule中
        :return:
        """
        self._molecule.add_bond(bond, self, other_atom)

    def link_bond_by_type(self, bond_type: str, other_atom):
        """ 通过molecule创建并连接bond，bond会被加入molecule

        :param bond_type: 待创建的键的类型
        :param other_atom: 本原子通过键连接的另一个原子，other_atom应该存在于molecule中
        :return:
        """
        self._molecule.add_bond_by_type(bond_type, self, other_atom)

    def remove_bond(self, bond) -> None:
        """ 移除键，并将移除过的键加入到已移除键的序列中

        :param bond: 待移除的键
        :return:
        """
        self._brock_bonds.append(bond)
        self._bonds.remove(bond)

    # ====================
    # -- Neighbor Atoms --

    def get_neighbors(self) -> Iterator:
        """ 获得临近原子

        :return: 临近原子的生成器
        """
        return self._molecule.get_neighbors(self)

    def get_neighbors_excepts(self, other_atoms: List) -> Iterator:
        """ 获得排除一些原子后的临近原子

        :param other_atoms: 需要排除的原子
        :return: 临近原子的生成器（除了需要排除的原子）
        """
        for atom in self.get_neighbors():
            if atom not in other_atoms:
                yield atom

    def get_neighbor_atoms_by_level(self, level: int) -> List:
        atoms = [self]
        except_atoms = [self]
        for i in range(level):
            next_level_atoms = []
            for atom in atoms:
                next_level_atoms.extend(list(atom.get_neighbors_excepts(except_atoms)))
            if i == level - 1:
                return next_level_atoms
            except_atoms.extend(next_level_atoms)
            atoms = next_level_atoms

    def is_neighbored(self, other) -> bool:
        for neighbor in self.get_neighbors():
            if neighbor == other:
                return True
        return False

    # ================
    # -- Properties --

    def set_property(self, key: str, value: object) -> None:
        self.properties[key] = value

    def get_property(self, key: str) -> Any:
        if key in self.properties.keys():
            return self.properties[key]
        else:
            return None

    def clear_property(self, key: str) -> None:
        if key in self.properties.keys():
            self.properties.pop(key)

    def is_dropped_in_rxn(self):
        return self.map_num is None or self.map_num == 0 or not self.is_visible

    def copy(self, save_parent=False, save_copy_hash=True):
        atom = Atom()
        self._copy_two_atom(atom, self, save_copy_hash)
        if save_parent:
            atom.set_parent(self)
        return atom

    def copy_atom(self, atom, without_map_num: bool = False):
        map_num = self.map_num
        self._copy_two_atom(self, atom, copy_hash=False)
        if without_map_num:
            self.map_num = map_num

    @staticmethod
    def _copy_two_atom(new_atom, old_atom, copy_hash: bool = True):
        new_atom.set_hybrid(old_atom.get_hybrid())
        new_atom.set_symbol(old_atom.get_symbol())
        new_atom.map_num = old_atom.map_num
        new_atom.is_visible = old_atom.is_visible
        new_atom.degree = old_atom.degree
        new_atom.properties = old_atom.properties.copy()
        new_atom.set_id(old_atom.get_id())
        new_atom.charge = old_atom.charge
        new_atom.is_aromatic = old_atom.is_aromatic
        new_atom.set_implicit_valence(old_atom.get_implicit_valence())
        new_atom.set_explicit_valence(old_atom.get_explicit_valence())
        new_atom.implicit_hs_num = old_atom.implicit_hs_num
        new_atom.explicit_hs_num = old_atom.explicit_hs_num
        if copy_hash:
            new_atom.set_copy_hash(old_atom.get_copy_hash())

    @staticmethod
    def instance_by_symbol(symbol: str):
        atom = Atom()
        atom.set_symbol(symbol)
        return atom
