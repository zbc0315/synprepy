# -*- coding: utf-8 -*-
# @Time     : 2020/4/20 16:22
# @Author   : zbc
# @Email    : zhangbc0315@outlook.com
# @FileName : molecule.py
# @Software : PyCharm

from typing import Iterator, List, Tuple, Dict, Any

import networkx as nx

from crkitpy.molgraph.atom import Atom
from crkitpy.molgraph.bond import Bond
from crkitpy.molgraph.atom_set import AtomSet
from crkitpy.draw.plt_draw import PltDraw
from crkitpy.const import ATOM_PROP
from crkitpy.const import BOND_TYPE
from crkitpy.chemdata.atomic_data import AtomicData


class Molecule:

    def __init__(self):
        self.mid = 0
        self.mol_graph = nx.Graph()
        self._mol_graph_for_match = None
        self._mol_graph_visible = None
        self._atoms: List[Atom] = []
        self._bonds: List[Bond] = []
        self._used_linked_bonds: List[Bond] = []
        self._used_broken_bonds: List[Bond] = []
        self._changed_charge: Dict[Atom, Tuple[int, int]] = {}
        self._changed_valence: Dict[Atom, Tuple[int, int]] = {}
        self._aromatic_cycles: List[AtomSet] = []
        self.properties = {}

    # def __del__(self):
    #     self._atoms.clear()
    #     self._bonds.clear()
    #     del self.mol_graph

    # ==========
    # -- copy --

    def copy(self, molecule=None, save_parent: bool = False, save_atom_copy_hash: bool = True, save_mid: bool = False) -> Tuple[Any, Dict[Atom, Atom], Dict[Bond, Bond]]:
        target_molecule = Molecule() if molecule is None else molecule
        atom_map = {}
        bond_map = {}
        # 复制原子
        for source_atom in self.get_atoms():
            atom_map[source_atom] = source_atom.copy(save_parent=save_parent, save_copy_hash=save_atom_copy_hash)
            if save_mid:
                atom_map[source_atom].set_property(ATOM_PROP.K_MOL_ID, self.mid)
            target_molecule.add_atom(atom_map[source_atom])
        # 复制键
        for source_bond in self.get_bonds():
            bond_map[source_bond] = source_bond.copy(save_parent=save_parent)
            target_molecule.add_bond(bond_map[source_bond],
                                     atom_map[source_bond.get_atoms()[0]],
                                     atom_map[source_bond.get_atoms()[1]])
        # 复制曾经断裂键的信息
        for source_bond in self.get_used_broken_bonds():
            target_molecule.add_used_broken_bond(bond_map[source_bond])
        # 复制曾经连接键的信息
        for source_v_bond in self.get_used_linked_bonds():
            target_v_bond = Bond.instance_by_bond_type(source_v_bond.get_bond_type())
            target_v_bond.add_atom(atom_map[source_v_bond.get_atoms()[0]])
            target_v_bond.add_atom(atom_map[source_v_bond.get_atoms()[1]])
            target_molecule.add_used_linked_bond(target_v_bond)
        # 复制电荷变化信息
        for source_atom, charges in self.get_changed_charge().items():
            target_molecule.add_changed_charge(atom_map[source_atom], charges[0], charges[1])
        # 复制价态变化信息
        for source_atom, valences in self.get_changed_valence().items():
            target_molecule.add_changed_valence(atom_map[source_atom], valences[0], valences[1])
        return target_molecule, atom_map, bond_map

    def copy_sub_mol(self, atoms: List[Atom]):
        sub_mol = Molecule()
        for atom in atoms:
            if atom not in self.get_atoms():
                raise ValueError('Cannot copy sub mol by atoms not belong this mol')
            sub_mol.add_atom(atom.copy(save_parent=True))
        bonds = []
        for i, begin_atom in enumerate(atoms[:-1]):
            for j, end_atom in enumerate(atoms[i+1:]):
                bond = self.get_bond(begin_atom, end_atom)
                if bond not in bonds and bond is not None:
                    bonds.append(bond)
        for bond in bonds:
            c_bond = bond.copy(save_parent=True)
            c_begin_atom = sub_mol.get_only_one_atom_by_copy_hash(bond.get_atoms()[0].get_copy_hash())
            c_end_atom = sub_mol.get_only_one_atom_by_copy_hash(bond.get_atoms()[1].get_copy_hash())
            sub_mol.add_bond(c_bond, c_begin_atom, c_end_atom)
        return sub_mol

    # ==================
    # -- add molecule --

    def add_molecule(self, mol, link_this_other_bonds: List[Tuple[Atom, str, Atom]], retain=False):
        if not retain:
            molecule, other_old_new_atoms, _ = mol.copy()
        else:
            molecule = mol
        molecule, this_old_new_atoms, _ = self.copy(molecule)
        for old_this_atom, bond_type, old_other_atom in link_this_other_bonds:
            new_this_atom = molecule.get_only_one_atom_by_copy_hash(old_this_atom.get_copy_hash())
            new_other_atom = molecule.get_only_one_atom_by_copy_hash(old_other_atom.get_copy_hash()) if not retain else old_other_atom
            # new_this_atom = this_old_new_atoms[old_this_atom]
            # new_other_atom = other_old_new_atoms[old_other_atom] if not retain else old_other_atom
            molecule.add_bond_by_type(bond_type, new_this_atom, new_other_atom)
        return molecule

    # =====================
    # -- union molecules --

    @classmethod
    def union_mols(cls, mols):
        molecule = Molecule()
        for i, mol in enumerate(mols):
            mol.mid = i
            molecule, _, __ = mol.copy(molecule=molecule, save_parent=True, save_mid=True)
        return molecule

    # ==========
    # -- draw --

    def draw(self, display=True):
        weights = [BOND_TYPE.BOND_WEIGHTS[bond.get_bond_type()] * 5 - 4 for bond in self.get_bonds_from_graph()]
        node_color = ['y' if atom.get_property(ATOM_PROP.K_NEED_SAVE) else 'w' for atom in self.get_atoms()]
        PltDraw.draw(self.mol_graph, width=weights, node_color=node_color)
        if display:
            PltDraw.show()

    # ====================
    # -- aromatic cycle --

    def add_aromatic_cycle(self, cycle: AtomSet) -> None:
        if cycle not in self._aromatic_cycles:
            self._aromatic_cycles.append(cycle)
            for atom in cycle.get_atoms():
                atom._add_aromatic_cycle(cycle)
                for bond in atom.get_bonds():
                    if bond.get_another_atom(atom) in cycle.get_atoms():
                        cycle.add_bond(bond)

    def get_aromatic_cycles(self) -> List[AtomSet]:
        return self._aromatic_cycles

    # ====================
    # -- changed charge --

    def add_changed_charge(self, atom: Atom, old_charge: int, new_charge: int):
        if atom not in self._changed_charge.keys():
            self._changed_charge[atom] = (old_charge, new_charge)

    def get_changed_charge(self) -> Dict[Atom, Tuple[int, int]]:
        return self._changed_charge

    # =====================
    # -- changed valence --

    def add_changed_valence(self, atom: Atom, old_valence: int, new_valence: int):
        if atom not in self._changed_valence.keys():
            self._changed_valence[atom] = (old_valence, new_valence)

    def get_changed_valence(self) -> Dict[Atom, Tuple[int, int]]:
        return self._changed_valence

    # ======================
    # -- used broken bond --

    def add_used_broken_bond(self, bond):
        self._used_broken_bonds.append(bond)

    def get_used_broken_bonds(self):
        return self._used_broken_bonds

    def clear_used_broken_bonds(self) -> None:
        self._used_broken_bonds = []

    # ======================
    # -- used linked bond --

    def add_used_linked_bond(self, bond):
        self._used_linked_bonds.append(bond)

    def add_used_linked_bond_by_type(self, bond_type, begin_atom, end_atom):
        bond = Bond.instance_by_bond_type(bond_type)
        bond.add_atom(begin_atom)
        bond.add_atom(end_atom)
        self.add_used_linked_bond(bond)

    def get_used_linked_bonds(self):
        return self._used_linked_bonds

    def clear_used_linked_bonds(self) -> None:
        self._used_linked_bonds = []

    # ==============
    # -- property --

    def set_property(self, key: str, value: Any) -> None:
        self.properties[key] = value

    def get_property(self, key: str) -> Any:
        if key in self.properties.keys():
            return self.properties[key]
        else:
            return None

    def clear_property(self, key: str) -> None:
        if key in self.properties.keys():
            self.properties.pop(key)

    def clear_atoms_property(self, key: str) -> None:
        for atom in self.get_atoms():
            atom.clear_property(key)

    # ==========
    # -- atom --

    def add_atom(self, atom: Atom) -> None:
        if not self.mol_graph.has_node(atom):
            self._atoms.append(atom)
            atom.set_molecule(self)
            self.mol_graph.add_node(atom, symbol=atom.get_symbol(), charge=atom.charge, valence=atom.get_valence(),
                                    is_aromatic=atom.is_aromatic,
                                    hs_num=atom.get_hs_num())

    def create_and_add_atom_by_symbol(self, symbol: str) -> Atom:
        atom = Atom.instance_by_symbol(symbol)
        self.add_atom(atom)
        return atom

    def get_atoms(self) -> List[Atom]:
        return self._atoms

    def get_atoms_num(self) -> int:
        return len(self.get_atoms())

    def get_atoms_by_prop(self, prop_key: str, prop_value: object) -> Iterator[Atom]:
        for atom in self.get_atoms():
            if atom.get_property(prop_key) == prop_value:
                yield atom

    def get_atom_by_id(self, aid: int) -> Atom:
        for atom in self.get_atoms():
            if atom.get_id() == aid:
                return atom

    def get_atoms_by_symbol(self, symbol: str) -> Iterator[Atom]:
        for atom in self.get_atoms():
            if atom.get_symbol() == symbol:
                yield atom

    def get_one_atom_by_map_num(self, map_num: int, filter_func = None) -> Atom:
        for atom in self.get_atoms_by_map_num(map_num):
            if filter_func is None or filter_func(atom):
                return atom

    def get_atoms_by_map_num(self, map_num: int) -> Iterator[Atom]:
        for atom in self.get_atoms():
            if atom.map_num == map_num:
                yield atom

    def get_atoms_by_map_nums(self, map_nums: List[int]) -> Iterator[Atom]:
        for atom in self.get_atoms():
            if atom.map_num in map_nums:
                yield atom

    def get_atoms_by_copy_hash(self, copy_hash: int) -> Iterator[Atom]:
        for atom in self.get_atoms():
            if atom.get_copy_hash() == copy_hash:
                yield atom

    def get_only_one_atom_by_copy_hash(self, copy_hash: int) -> Atom:
        atoms = list(self.get_atoms_by_copy_hash(copy_hash))
        if len(atoms) != 1:
            raise ValueError('Must contain one atom have this copy_hash, but get {}'.format(len(atoms)))
        return atoms[0]

    def get_atom_by_parent(self, parent_atom) -> Atom:
        for atom in self._atoms:
            if atom.get_parent() == parent_atom:
                return atom

    def get_neighbors(self, atom: Atom) -> Iterator[Atom]:
        # TODO 检查atom是否在mol_graph中
        return self.mol_graph.neighbors(atom)

    def get_atoms_by_specific_property(self, prop_key, prop_value) -> Iterator[Atom]:
        for atom in self.get_atoms():
            if atom.get_property(prop_key) == prop_value:
                yield atom

    def get_reacted_atom(self) -> Iterator[Atom]:
        return self.get_atoms_by_specific_property(ATOM_PROP.K_IS_REACTED, True)

    def get_saved_atom(self) -> Iterator[Atom]:
        return self.get_atoms_by_specific_property(ATOM_PROP.K_NEED_SAVE, True)

    # ==========
    # -- bond --

    def add_bond(self, bond: Bond, begin_atom: Atom, end_atom: Atom) -> None:
        if bond.bond_to_r==False and (begin_atom.get_symbol() == '*' or end_atom.get_symbol() == '*'):
            bond.bond_to_r = True
        begin_atom.add_bond(bond)
        end_atom.add_bond(bond)
        self._bonds.append(bond)
        bond.set_molecule(self)
        self.mol_graph.add_edge(begin_atom, end_atom, bond_type=bond.get_bond_type(), bond_to_r=bond.bond_to_r,
                                bond_weight=BOND_TYPE.BOND_WEIGHTS[bond.get_bond_type()])
        # 在molhandler/cycle_detector中也设置了键的参数

    def add_bond_by_type(self, bond_type: str, begin_atom: Atom, end_atom: Atom, bond_to_r: bool = False) -> None:
        # TODO 检查begin_atom与end_atom之间是否已经有bond
        bond = Bond.instance_by_bond_type(bond_type)
        if begin_atom.get_symbol() == '*' or end_atom.get_symbol() == '*':
            bond.bond_to_r = True
        self.add_bond(bond, begin_atom, end_atom)

    def get_bond(self, begin_atom: Atom, end_atom: Atom) -> Bond:
        for bond in begin_atom.get_bonds():
            for atom in bond.get_atoms():
                if atom == end_atom:
                    return bond

    def get_bonds(self) -> List[Bond]:
        bonds = []
        for atom in self.get_atoms():
            for b in atom.get_bonds():
                if b not in bonds:
                    bonds.append(b)
        return bonds

    def get_bonds_from_mol(self) -> List[Bond]:
        return self._bonds

    def get_bonds_from_graph(self) -> List[Bond]:
        atom_pairs: List[Tuple[Atom, Atom]] = self.mol_graph.edges
        return [atom_pair[0].get_bond(atom_pair[1]) for atom_pair in atom_pairs]

    def get_bonds_by_copy_hash(self, copy_hash: int) -> Iterator[Bond]:
        for bond in self.get_bonds():
            if bond.get_copy_hash() == copy_hash:
                yield bond

    def get_only_one_bond_by_copy_hash(self, copy_hash: int) -> Bond:
        bonds = list(self.get_bonds_by_copy_hash(copy_hash))
        if len(bonds) != 1:
            raise ValueError('Must contain one bond have this copy hash')
        return bonds[0]

    def remove_bond(self, bond: Bond) -> None:
        self.mol_graph.remove_edge(bond.get_atoms()[0], bond.get_atoms()[1])
        self._bonds.remove(bond)
        bond.get_atoms()[0].remove_bond(bond)
        bond.get_atoms()[1].remove_bond(bond)

    def remove_atom(self, atom: Atom) -> None:
        neighbor_atoms = list(self.get_neighbors(atom))
        for neighbor_atom in neighbor_atoms:
            bond = self.get_bond(atom, neighbor_atom)
            neighbor_atom.remove_bond(bond)
        self.mol_graph.remove_node(atom)
        self._atoms.remove(atom)

    def remove_group_by_copy_hash(self, group):
        for atom in group.get_atoms():
            self.remove_atom(self.get_only_one_atom_by_copy_hash(atom.get_copy_hash()))

    def get_mol_graph_for_match(self):
        if self._mol_graph_for_match is None:
            self._mol_graph_for_match = nx.Graph()
            for atom in self.get_atoms():
                if atom.get_symbol() == '*':
                    continue
                self._graph_add_atom_for_template(self._mol_graph_for_match, atom)
            for bond in self.get_bonds():
                if bond.is_link_atom_by_symbol('*'):
                    continue
                self._graph_add_bond_for_template(self._mol_graph_for_match, bond)
        return self._mol_graph_for_match

    # =======================
    # -- mol graph visible --

    def get_mol_graph_visible(self):
        if self._mol_graph_visible is None:
            self._mol_graph_visible = nx.Graph()
            for atom in self.get_atoms():
                if atom.is_visible:
                    self._graph_add_atom_for_template(self._mol_graph_visible, atom)
                for bond in self.get_bonds():
                    if all([a.is_visible for a in bond.get_atoms()]):
                        self._graph_add_bond_for_template(self._mol_graph_visible, bond)
        return self._mol_graph_visible

    @classmethod
    def _graph_add_atom_for_template(cls, graph: nx.Graph, atom: Atom) -> None:
        require_valence = atom.get_property(ATOM_PROP.K_REQUIRE_VALENCE)
        require_valence = require_valence if require_valence is not None else 0

        max_valence = atom.get_property(ATOM_PROP.K_MAX_VALENCE)
        if max_valence is None:
            max_valence = max(AtomicData.get_valences(atom.get_symbol()))
        graph.add_node(atom, symbol=atom.get_symbol(), valence=atom.get_valence(), charge=atom.charge,
                       is_aromatic=atom.is_aromatic, no_h_bond_valence=atom.get_bonds_weight_without_h(),
                       require_valence=require_valence,
                       max_valence=max_valence,
                       hs_num=atom.get_hs_num())

    @classmethod
    def _graph_add_bond_for_template(cls, graph: nx.Graph, bond: Bond) -> None:
        graph.add_edge(bond.get_atoms()[0], bond.get_atoms()[1], bond_type=bond.get_bond_type(), bond_to_r=bond.bond_to_r)

    def clear_mapping(self):
        for atom in self.get_atoms():
            atom.map_num = None

    def compose_graph(self, other_mol) -> None:
        self.mol_graph = nx.compose(self.mol_graph, other_mol.mol_graph)
        for atom in other_mol.get_atoms():
            self._atoms.append(atom)
        for bond in other_mol.get_bonds():
            self._bonds.append(bond)

    def get_detached_sub_mols(self) -> Iterator:
        for atoms in nx.connected_components(self.mol_graph):
            yield self.copy_sub_mol(list(atoms))

    def detect_aromatic_ring(self):
        di_graph = nx.DiGraph(self.mol_graph)
        cycles = list(nx.simple_cycles(di_graph))
        i = 1
