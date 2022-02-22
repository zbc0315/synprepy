# -*- coding: utf-8 -*-
# @Time     : 2020/5/29 9:39
# @Author   : zbc
# @Email    : zhangbc0315@outlook.com
# @FileName : rxn_reverse_template_extractor.py
# @Software : PyCharm

from typing import List, Iterator, Tuple, Sequence

from crkitpy.rxnhandler.rxn_center_detector import RxnCenterDetector
from crkitpy.rxnhandler.rxn_template_extractor import RxnTemplateExtractor
from crkitpy.molgraph.reaction import Reaction
from crkitpy.molgraph.molecule import Molecule, AtomSet
from crkitpy.molgraph.atom import Atom
from crkitpy.molgraph.bond import Bond
from crkitpy.const import ATOM_PROP
from crkitpy.molhandler.aromatic_cycle_detector import AromaticCycleDetector
from crkitpy.smilesparse.smiles_writer import SmilesWriter


class RxnReverseTemplateExtractor:
    """
    map num 在反应物与产物中皆应只出现一次

    """

    @classmethod
    def extract_rxn_templates(cls, reaction: Reaction,
                              highest_react_level: int = 1,
                              with_aromatic: bool = True) -> Iterator[Molecule]:
        """
        :param reaction:
        :param highest_react_level:
        :param with_aromatic:
        :return:
        """
        smis = []
        for molecule in reaction.get_reactants():
            molecule.draw()
            smis.append(SmilesWriter.mol_to_smiles(molecule, need_reset_map_num=False, need_normal_smiles=False))
        for molecule in reaction.get_products():
            molecule.draw()
            smis.append(SmilesWriter.mol_to_smiles(molecule, need_reset_map_num=False, need_normal_smiles=False))

        for molecule in reaction.get_molecules():
            molecule.clear_atoms_property(ATOM_PROP.K_REACT_LEVEL)
        RxnTemplateExtractor.check_available_rxn(reaction)
        if with_aromatic:
            RxnTemplateExtractor.detect_aromatic_cycles(reaction)
        RxnCenterDetector.detect_rxn_center(reaction)
        for react_level in range(1, highest_react_level + 1):
            RxnTemplateExtractor.sign_reactant_react_level_in_rxn(reaction, react_level, with_aromatic)
        # RxnTemplateExtractor.check_one_step_rxn_by_react_level_and_prop(reaction)  # TODO 目前不考虑产物中有多个重复 map num的情况
        for rxn_template in cls._extract_core_template(reaction):
            cls._add_dropped_atoms(reaction, rxn_template)
            cls._detect_used_linked_bonds(reaction, rxn_template)
            cls._detect_used_broken_bonds_and_changed_props(reaction, rxn_template)
            yield rxn_template

    # ==========================
    # -- detect changed props --

    @classmethod
    def _detect_changed_props(cls, rxn_template: Molecule, p_atom: Atom, r_atom: Atom) -> None:
        if p_atom.get_valence() != r_atom.get_valence():
            rxn_template.add_changed_valence(p_atom, r_atom.get_valence(), p_atom.get_valence())
        if p_atom.charge != r_atom.charge:
            rxn_template.add_changed_charge(p_atom, r_atom.charge, p_atom.charge)

    # ==============================
    # -- detect used broken bonds --

    @classmethod
    def _detect_used_broken_bonds_and_changed_props(cls, reaction: Reaction, rxn_template: Molecule):
        for p_atom in rxn_template.get_atoms_by_prop(ATOM_PROP.K_REACT_LEVEL, 0):
            r_atom = cls._get_mapping_atom_from_reactants(reaction, p_atom.map_num)
            cls._detect_changed_props(rxn_template, p_atom, r_atom)
            for p_n_atom in p_atom.get_neighbors():
                p_bond = p_atom.get_bond(p_n_atom)
                r_n_atom = cls._get_mapping_atom_from_reactants(reaction, p_n_atom.map_num)
                if r_n_atom is None:
                    raise ValueError('Cannot get mapping atom from reactant by product atom')
                r_bond = r_atom.get_bond(r_n_atom)
                if r_bond is None or r_bond.get_bond_type() != p_bond.get_bond_type():
                    if p_bond not in rxn_template.get_used_broken_bonds():
                        rxn_template.add_used_broken_bond(p_bond)

    # ==============================
    # -- detect used linked bonds --

    @classmethod
    def _detect_used_linked_bonds(cls, reaction: Reaction, rxn_template: Molecule) -> None:

        for p_atom in rxn_template.get_atoms_by_prop(ATOM_PROP.K_REACT_LEVEL, 0):
            r_atom = cls._get_mapping_atom_from_reactants(reaction, p_atom.map_num)
            for r_n_atom in r_atom.get_neighbors():
                if r_n_atom.is_dropped_in_rxn():
                    rxn_template.add_used_linked_bond_by_type(r_atom.get_bond(r_n_atom).get_bond_type(),
                                                              p_atom, r_n_atom)
                else:
                    p_n_atom = cls._get_mapping_atom_from_atoms(p_atom.get_neighbors(), r_n_atom.map_num)
                    r_bond = r_atom.get_bond(r_n_atom)
                    if p_n_atom is None:
                        p_map_n_atom = cls._get_mapping_atom_from_atoms(rxn_template.get_atoms(), r_n_atom.map_num)
                        rxn_template.add_used_linked_bond_by_type(r_bond.get_bond_type(), p_atom, p_map_n_atom)
                    else:
                        p_bond = p_atom.get_bond(p_n_atom)
                        if p_bond.get_bond_type() != r_bond.get_bond_type():
                            rxn_template.add_used_linked_bond_by_type(r_bond.get_bond_type(), p_atom, p_n_atom)
        cls._duplicate_and_check_used_linked_bonds(rxn_template)

    @classmethod
    def _duplicate_and_check_used_linked_bonds(cls, rxn_template: Molecule) -> None:
        bonds = rxn_template.get_used_linked_bonds().copy()
        rxn_template.clear_used_linked_bonds()
        added_bonds = []
        for bond in bonds:
            if bond in added_bonds:
                continue
            if any([atom.is_dropped_in_rxn() for atom in bond.get_atoms()]):
                rxn_template.add_used_linked_bond(bond)
            else:
                equal_bonds = cls._get_equal_bonds_by_type_and_atoms(bonds, bond)
                if len(equal_bonds) == 2:
                    added_bonds.extend(equal_bonds)
                    rxn_template.add_used_linked_bond(bond)
                elif len(equal_bonds) > 2:
                    raise ValueError('something unknown wrong')

    @staticmethod
    def _get_equal_bonds_by_type_and_atoms(bonds: List[Bond], bond: Bond) -> List[Bond]:
        equal_bonds = []
        for _bond in bonds:
            if bond.equal_in_type_and_atoms(_bond):
                equal_bonds.append(_bond)
        return equal_bonds

    @classmethod
    def _get_mapping_atom_from_atoms(cls, atoms: Iterator[Atom], map_num: int) -> Atom:
        for atom in atoms:
            if atom.map_num == map_num:
                return atom

    @classmethod
    def _get_mapping_atom_from_reactants(cls, reaction: Reaction, map_num: int) -> Atom:
        for reactant in reaction.get_reactants():
            for atom in reactant.get_atoms_by_map_num(map_num):
                return atom

    # =======================
    # -- add dropped atoms --

    @classmethod
    def _add_dropped_atoms(cls, reaction: Reaction, rxn_template: Molecule):
        for reactant in reaction.get_reactants():
            invisible_atoms = []
            for atom in reactant.get_atoms():
                if atom.is_dropped_in_rxn():
                    atom.is_visible = False
                    invisible_atoms.append(atom)
                    rxn_template.add_atom(atom)
            for bond in reactant.get_bonds():
                if all([a in invisible_atoms for a in bond.get_atoms()]): # TODO 可以简化
                    rxn_template.add_bond_by_type(bond.get_bond_type(), bond.get_atoms()[0], bond.get_atoms()[1])

    # ======================
    # -- extract rxn area --

    @classmethod
    def _extract_core_template(cls, reaction: Reaction) -> Iterator[Molecule]:
        for product in reaction.get_products():
            atoms = []
            for atom in product.get_atoms():
                if atom.get_property(ATOM_PROP.K_REACT_LEVEL) is not None:
                    atoms.append(atom)
            yield product.copy_sub_mol(atoms)

    # ==========================================
    # ------------------ DROP ------------------

    @classmethod
    def _expand_atom(cls, atom: Atom) -> None:
        """ 扩展反应区域
        1. 本原子在芳环中：芳环上其他反应等级为None的原子，设置其反应等级为本原子的下一级
        2. 本原子临近的原子在芳环上：芳环上反应等级为None的原子，设置其反应等级为本原子的下一级
        3. 本原子临近的所有反应等级为None的原子，设置其反应等级为本原子的下一级

        :param atom:
        :return:
        """
        level = atom.get_property(ATOM_PROP.K_REACT_LEVEL)
        next_level = level + 1
        if atom.is_aromatic:
            for cycle in atom.get_aromatic_cycles():
                cls._set_aromatic_cycle_react_level(cycle, next_level)
        for neighbor in atom.get_neighbors():
            if neighbor.get_property(ATOM_PROP.K_REACT_LEVEL) is None:
                if neighbor.is_aromatic:
                    for cycle in neighbor.get_aromatic_cycles():
                        cls._set_aromatic_cycle_react_level(cycle, next_level)
                else:
                    neighbor.set_property(ATOM_PROP.K_REACT_LEVEL, next_level)

    # ====================================
    # -- set aromatic cycle react level --

    @classmethod
    def _set_aromatic_cycle_react_level(cls, cycle: AtomSet, react_level: int) -> None:
        """ 将芳环上所有react level为None的都设置为react_level

        :param cycle:
        :param react_level:
        :return:
        """
        for atom in cycle.get_atoms():
            if atom.get_property(ATOM_PROP.K_REACT_LEVEL) is None:
                atom.set_property(ATOM_PROP.K_REACT_LEVEL, react_level)

    # ==========================
    # -- unity aromatic cycle --

    @classmethod
    def _unity_aromatic_in_rxn(cls, reaction: Reaction, unity_level: int) -> None:
        """ 使芳环上的原子团结一致
        芳环上有一个原子的反应等级等于react_level，则该芳环上上所有原子的反应等级都是react_level
        注意：认为待处理芳环只有两种可能：1. 芳环上所有原子的反应等级都相同
                                    2. 芳环上部分原子的反应等级为None，部分原子的反应等级都等于react_level
             请避免出现第三种情况
        :param reaction:
        :param unity_level:
        :return:
        """
        for product in reaction.get_products():
            cls._unity_aromatic_in_mol(product, unity_level)

    @classmethod
    def _unity_aromatic_in_mol(cls, molecule: Molecule, unity_level: int) -> None:
        """ 使芳环上的原子团结一致
        芳环上有一个原子的反应等级等于react_level，则该芳环上上所有原子的反应等级都是react_level
        注意：认为待处理芳环只有两种可能：1. 芳环上所有原子的反应等级都相同
                                    2. 芳环上部分原子的反应等级为None，部分原子的反应等级都等于react_level
             请避免出现第三种情况
        :param molecule:
        :param unity_level:
        :return:
        """
        parsed_cycles = []
        for atom in molecule.get_atoms_by_prop(ATOM_PROP.K_IS_REACTED, True):
            for cycle in atom.get_aromatic_cycles():
                if cycle not in parsed_cycles:
                    cls._unity_aromatic_cycle(cycle, unity_level)
                    parsed_cycles.append(cycle)

    @staticmethod
    def _unity_aromatic_cycle(cycle: AtomSet, unity_level: int) -> None:
        """ 使芳环上的原子团结一致
        芳环上有一个原子的反应等级等于react_level，则该芳环上上所有原子的反应等级都是react_level
        注意：认为待处理芳环只有两种可能：1. 芳环上所有原子的反应等级都相同
                                    2. 芳环上部分原子的反应等级为None，部分原子的反应等级都等于react_level
             请避免出现第三种情况
        :param molecule:
        :param unity_level:
        :return:
        """
        all_other_level = True
        other_level = None
        need_unity = False
        # 检查
        for atom in cycle.get_atoms():
            level = atom.get_property(ATOM_PROP.K_REACT_LEVEL)
            if level is None:
                all_other_level = False
            elif level != unity_level:
                if not all_other_level:
                    # 芳环上出现None与其他level，错误
                    raise ValueError("aromatic unity error")
                if other_level is None:
                    other_level = level
                elif other_level != level:
                    # 芳环上出现不同的level，错误
                    raise ValueError("aromatic unity error")
            elif other_level is not None:
                # 芳环上出现不同的level，错误
                raise ValueError("aromatic unity error")
            else:
                need_unity = True

        # 赋值
        if need_unity:
            for atom in cycle.get_atoms():
                atom.set_property(ATOM_PROP.K_REACT_LEVEL, unity_level)

