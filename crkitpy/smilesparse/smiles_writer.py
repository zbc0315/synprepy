# -*- coding: utf-8 -*-
# @Time     : 2020/4/21 14:24
# @Author   : zbc
# @Email    : zhangbc0315@outlook.com
# @FileName : smiles_writer.py
# @Software : PyCharm

from typing import List, Dict

from rdkit.Chem import AllChem

from crkitpy.molgraph.molecule import Molecule, Atom
from crkitpy.const import ATOM_PROP, BOND_TYPE
from crkitpy.rxnhandler.map_num_reset import MapNumReset
from crkitpy.molhandler.molecule_checker import MoleculeChecker
from crkitpy.molhandler.aromatic_handler import AromaticHandler
from crkitpy.molhandler.molecule_hs_num_corrector import MoleculeHsNumCorrector
from crkitpy.chemdata.atomic_data import AtomicData
from crkitpy.errors.smiles_error import SmilesWriterError


class SmilesWriter:
    """ Smiles是本程序定义的，类似于smiles，是template的唯一字符串表示
        TODO 未处理分子中的氢
    """

    def __init__(self):
        self.first_atom = None
        self.passed_atoms: List[Atom] = []
        self.parent_atoms: Dict[Atom, Atom] = {}
        self.childs_atoms: Dict[Atom, List[Atom]] = {}
        self.cycle_num = 0
        self.branch_start_atoms = []
        self.parsed_cycle_nums = []
        self.need_normal_smiles = True

    def _atom_to_smiles(self, atom: Atom) -> str:
        need_bracket = self._need_bracket(atom)
        symbol = atom.get_symbol()
        if atom.is_aromatic and not self.need_normal_smiles:
            symbol = symbol.lower()
        atom_smiles_list = [symbol]
        if atom.get_hs_num() > 0 and need_bracket:
            atom_smiles_list.append(self._get_hs_symbol(atom.get_hs_num()))
        if atom.charge is not None and atom.charge != 0:
            atom_smiles_list.append(self._get_charge_symbol(atom.charge))
        if atom.map_num is not None and atom.map_num > 0:
            atom_smiles_list.append(self._get_map_num_symbol(atom.map_num))
        if need_bracket:
            atom_smiles_list.insert(0, '[')
            atom_smiles_list.append(']')
        if atom.get_property(ATOM_PROP.K_CYCLE_NUM) is not None:
            atom_smiles_list.append(self._get_cycle_nums_symbol(atom.get_property(ATOM_PROP.K_CYCLE_NUM)))
        return ''.join(atom_smiles_list)

    def _need_bracket(self, atom: Atom) -> bool:
        if self.need_normal_smiles:
            valences = AtomicData.get_valences(atom.get_symbol())
            return len(atom.get_symbol()) > 1 \
                   or atom.charge != 0 \
                   or (atom.map_num is not None and atom.map_num != 0) \
                   or atom.get_valence() not in valences
        else:
            return True

    @staticmethod
    def _get_hs_symbol(hs_num: int) -> str:
        if hs_num == 1:
            return 'H'
        elif hs_num > 1:
            return f'H{hs_num}'

    @staticmethod
    def _get_charge_symbol(charge: int) -> str:
        if charge == -1:
            return '-'
        elif charge == +1:
            return '+'
        elif charge < 0:
            return str(charge)
        elif charge > 0:
            return f'+{charge}'

    @staticmethod
    def _get_map_num_symbol(map_num: int) -> str:
        return f':{map_num}'

    def _get_cycle_nums_symbol(self, cycle_nums: Dict[int, str]) -> str:
        cycle_num_symbol_sb = []
        sep = ''
        if any([num > 9 for num in cycle_nums.keys()]):
            sep = '%'
        for num, bond_type in cycle_nums.items():
            if bond_type in [BOND_TYPE.SINGLE, BOND_TYPE.AROMATIC] or num in self.parsed_cycle_nums:
                cycle_num_symbol_sb.append(str(num))
            else:
                bond_symbol = BOND_TYPE.BOND_TO_SYMBOLS[bond_type]
                cycle_num_symbol_sb.append(f'{bond_symbol}{num}')
                self.parsed_cycle_nums.append(num)
        return sep.join(cycle_num_symbol_sb)

    # ========================
    # -- template to smiles --

    def _mol_to_smiles(self, molecule: Molecule,
                       need_reset_map_num: bool = True) -> str:
        if self.need_normal_smiles:
            AromaticHandler.remove_aromatic(molecule)
        molecule.clear_atoms_property(ATOM_PROP.K_CYCLE_NUM)
        if need_reset_map_num:
            MapNumReset.reset_molecule_map_num(molecule)
        smileses = []
        for sub_mol in molecule.get_detached_sub_mols():
            smileses.append(self._sub_mol_to_smiles(sub_mol))
        smileses.sort()
        return '.'.join(smileses)

    def _sub_mol_to_smiles(self, molecule: Molecule) -> str:
        self.first_atom = self._get_first_atom(molecule)
        self.branch_start_atoms = [self.first_atom]
        self.passed_atoms = [self.first_atom]
        mol_smiles_list = [self.first_atom]
        bracket_num = 0
        bond_num = 0
        while True:
            branch_start_atoms = self.branch_start_atoms.copy()
            for branch_start_atom in branch_start_atoms:
                branch_idx = mol_smiles_list.index(branch_start_atom)
                exclude_atoms = self._get_passed_neighbors(branch_start_atom)
                if self._has_no_branch(branch_start_atom, exclude_atoms):
                    self.branch_start_atoms.remove(branch_start_atom)
                    break
                branch = self._get_branch(branch_start_atom, exclude_atoms)
                if len(branch) > 0:
                    if branch_idx == len(mol_smiles_list) - 1:
                        mol_smiles_list.extend(branch)
                    else:
                        bracket_num += 2
                        mol_smiles_list = mol_smiles_list[:branch_idx + 1] + ['('] + branch + [')'] + mol_smiles_list[
                                                                                                      branch_idx + 1:]
                        branch_idx += 1
                    bond = branch_start_atom.get_bond(branch[0])
                    if bond.get_bond_type() not in [BOND_TYPE.SINGLE, BOND_TYPE.AROMATIC]:
                        branch_idx += 1
                        mol_smiles_list.insert(branch_idx, bond.symbol)  # TODO
                        bond_num += 1
            if len(self.branch_start_atoms) == 0:
                break
        mol_smiles_str_list = []
        for a in mol_smiles_list:
            if isinstance(a, Atom):
                mol_smiles_str_list.append(self._atom_to_smiles(a))
            else:
                mol_smiles_str_list.append(a)
        return ''.join(mol_smiles_str_list)

    @staticmethod
    def _get_first_atom(molecule: Molecule) -> Atom:
        min_id = None
        min_id_atom = None
        for atom in molecule.get_atoms():
            if atom.get_id() is None:
                molecule.draw()
                print(1)
            if min_id is None or min_id > atom.get_id():
                min_id_atom = atom
                min_id = atom.get_id()
        return min_id_atom

    def _has_no_branch(self, atom: Atom, exclude_atoms: List[Atom]) -> bool:
        neighbor_num = len(list(atom.get_neighbors()))
        if neighbor_num == len(exclude_atoms):
            return True
        if atom.get_property(ATOM_PROP.K_CYCLE_NUM) is not None and \
                neighbor_num == len(exclude_atoms) + len(atom.get_property(ATOM_PROP.K_CYCLE_NUM)):
            return True
        return False

    def _get_branch(self, start_atom: Atom, exclude_atoms: List[Atom] = None) -> List[Atom]:
        branch: List[Atom] = []
        atom = start_atom
        while True:
            sorted_next_atoms = self._get_sorted_next_atoms(atom, exclude_atoms)
            next_atom = None
            for n_atom in sorted_next_atoms:
                if n_atom in self.passed_atoms:
                    self._sign_cycle_num(atom, n_atom)
                elif next_atom is not None:
                    if atom not in self.branch_start_atoms:
                        self.branch_start_atoms.append(atom)
                    break
                else:
                    next_atom = n_atom
            if next_atom is None:
                break
            else:
                bond = atom.get_bond(next_atom)
                if len(branch) > 0 and bond.get_bond_type() not in [BOND_TYPE.SINGLE, BOND_TYPE.AROMATIC]:
                    branch.append(bond.symbol)
                branch.append(next_atom)
                self.parent_atoms[next_atom] = atom
                self._add_child(atom, next_atom)
                self.passed_atoms.append(next_atom)
            atom = next_atom
            exclude_atoms = None
        return branch

    def _sign_cycle_num(self, begin_atom: Atom, end_atom: Atom):
        self.cycle_num += 1
        begin_cycle_nums = begin_atom.get_property(ATOM_PROP.K_CYCLE_NUM)
        end_cycle_nums = end_atom.get_property(ATOM_PROP.K_CYCLE_NUM)
        if begin_cycle_nums is not None and end_cycle_nums is not None:
            inter_cycle_nums = set(begin_cycle_nums.keys()).intersection(set(end_cycle_nums.keys()))
            if len(inter_cycle_nums) != 0:
                return
        bond_type = begin_atom.get_bond(end_atom).get_bond_type()
        self._add_cycle_num(begin_atom, self.cycle_num, bond_type)
        self._add_cycle_num(end_atom, self.cycle_num, bond_type)

    def _add_cycle_num(self, atom: Atom, num: int, bond_type: str) -> None:
        if atom.get_property(ATOM_PROP.K_CYCLE_NUM) is None:
            atom.set_property(ATOM_PROP.K_CYCLE_NUM, {})
        atom.get_property(ATOM_PROP.K_CYCLE_NUM)[num] = bond_type

    def _get_sorted_next_atoms(self, atom: Atom, exclude_atoms: List[Atom] = None) -> List[Atom]:
        nbr_atoms = list(atom.get_neighbors())
        if len(self.parent_atoms) > 0 and atom != self.first_atom:
            nbr_atoms.remove(self.parent_atoms[atom])
        if exclude_atoms is not None:
            for atom in exclude_atoms:
                if atom in nbr_atoms:
                    nbr_atoms.remove(atom)
        nbr_atoms.sort(key=lambda x: x.get_id() if x.get_id() is not None else x.map_num)
        return nbr_atoms

    def _add_child(self, atom: Atom, child_atom: Atom) -> None:
        if atom not in self.childs_atoms.keys():
            self.childs_atoms[atom] = []
        self.childs_atoms[atom].append(child_atom)

    def _get_passed_neighbors(self, atom) -> List[Atom]:
        exclude_atoms = []
        if atom in self.parent_atoms:
            exclude_atoms.append(self.parent_atoms[atom])
        if atom in self.childs_atoms:
            exclude_atoms.extend(self.childs_atoms[atom])
        return exclude_atoms

    @classmethod
    def mol_to_smiles(cls, molecule: Molecule,
                      need_reset_map_num: bool = True,
                      need_normal_smiles: bool = False,
                      clear_map_num: bool = None) -> str:
        if len(molecule.get_atoms()) == 0:
            raise SmilesWriterError('wrong molecule with 0 atoms')
        if clear_map_num is None:
            clear_map_num = need_normal_smiles
        if clear_map_num:
            molecule.clear_mapping()
        # if need_normal_smiles:
        #     MoleculeChecker.check_molecule(molecule, need_correct_hs_num=True)
        # else:
        #     MoleculeHsNumCorrector.correct_hs_num(molecule)
        smiles_writer = SmilesWriter()
        smiles_writer.need_normal_smiles = need_normal_smiles
        return smiles_writer._mol_to_smiles(molecule, need_reset_map_num=need_reset_map_num)

    # =====================
    # -- inchi to smiles --

    @classmethod
    def inchi_to_smiles(cls, inchi: str) -> str:
        rdmol = AllChem.MolFromInchi(inchi)
        return AllChem.MolToSmiles(rdmol)
