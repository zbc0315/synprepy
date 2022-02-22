# -*- coding: utf-8 -*-
# @Time     : 2020/6/4 14:10
# @Author   : zbc
# @Email    : zhangbc0315@outlook.com
# @FileName : rxn_template_extractor.py
# @Software : PyCharm

from typing import Dict, List, Tuple

from crkitpy.molgraph.reaction import Reaction, Molecule
from crkitpy.molgraph.atom import Atom
from crkitpy.molhandler.aromatic_cycle_detector import AromaticCycleDetector
from crkitpy.rxnhandler.rxn_center_detector import RxnCenterDetector
from crkitpy.const import ATOM_PROP
from crkitpy.errors.rxn_error import *


class RxnTemplateExtractor:

    @classmethod
    def extract_rxn_template(cls, reaction: Reaction,
                             highest_react_level: int = 1,
                             with_aromatic: bool = True) -> Reaction:
        """ 抽取反应模板

        :param reaction: 化学反应
        :param highest_react_level: 在反应模板中需要考虑此反应等级及以下的原子，0级原子为参与反应的原子，1级原子与0级原子相邻，但不参与反应，以此类推
        :param with_aromatic: 是否需要统一设置芳环上的原子的反应等级
        :return:
        """
        cls.check_available_rxn(reaction)
        RxnCenterDetector.detect_rxn_center(reaction)
        if with_aromatic:
            cls.detect_aromatic_cycles(reaction)
            cls._sign_reactant_aromatic_cycle(reaction)
        for react_level in range(1, highest_react_level + 1):
            cls.sign_reactant_react_level_in_rxn(reaction, react_level, with_aromatic)
        # cls.check_one_step_rxn_by_react_level_and_prop(reaction) # TODO 在check_available_rxn中已经检查要求产物中不得出现重复map_num
        cls._sign_product_react_level(reaction)
        return cls._union_reactants_and_products(cls._extract_template_from_rxn(reaction))

    # ==============================
    # -- sign product react level --

    @classmethod
    def _sign_product_react_level(cls, rxn: Reaction) -> None:
        for reactant in rxn.get_reactants():
            for atom in reactant.get_atoms():
                if atom.get_property(ATOM_PROP.K_REACT_LEVEL) is None:
                    continue
                p_atom = cls._get_one_mapped_atom_from_mols(rxn.get_products(), atom.map_num)
                if p_atom is None:
                    continue
                p_atom.set_property(ATOM_PROP.K_REACT_LEVEL, atom.get_property(ATOM_PROP.K_REACT_LEVEL))

    # =========================
    # -- sign aromatic cycle --

    @classmethod
    def _sign_reactant_aromatic_cycle(cls, rxn: Reaction):
        for mol in rxn.get_reactants():
            for atom in mol.get_atoms_by_prop(ATOM_PROP.K_REACT_LEVEL, 0):
                cls._sign_react_level_in_aromatic_cycles(atom, 1)

    # ====================================
    # -- get mapped atom from molecules --

    @classmethod
    def _get_one_mapped_atom_from_mols(cls, mols: List[Molecule], map_num: int) -> Atom:
        for mol in mols:
            for atom in mol.get_atoms():
                if atom.map_num == map_num:
                    return atom

    # ==================================
    # -- union reactants and products --

    @classmethod
    def _union_reactants_and_products(cls, reaction: Reaction) -> Reaction:
        new_reaction = Reaction()
        new_reaction.add_reactant(Molecule.union_mols(reaction.get_reactants()))
        new_reaction.add_product(Molecule.union_mols(reaction.get_products()))
        return new_reaction

    # ======================
    # -- extract template --

    @classmethod
    def _extract_template_from_rxn(cls, reaction: Reaction) -> Reaction:
        template = Reaction()
        for reactant in reaction.get_reactants():
            template.add_reactant(cls._extract_template_from_mol(reactant))
        for product in reaction.get_products():
            template.add_product(cls._extract_template_from_mol(product))
        return template

    @classmethod
    def _extract_template_from_mol(cls, molecule: Molecule) -> Molecule:
        atoms = []
        for atom in molecule.get_atoms():
            if atom.get_property(ATOM_PROP.K_REACT_LEVEL) is not None or atom.is_dropped_in_rxn():
                atoms.append(atom)
        return molecule.copy_sub_mol(atoms)

    # ======================
    # -- sign react level --

    @classmethod
    def sign_reactant_react_level_in_rxn(cls, reaction: Reaction, react_level: int, with_aromatic: bool) -> None:
        for molecule in reaction.get_reactants():
            cls._sign_react_level_in_mol(molecule, react_level, with_aromatic)

    @classmethod
    def _sign_react_level_in_mol(cls, molecule: Molecule, react_level: int, with_aromatic: bool) -> None:
        if react_level < 0:
            raise ValueError(f'can not sign react_level: {react_level}')
        lower_react_level = react_level - 1
        for atom in molecule.get_atoms_by_prop(ATOM_PROP.K_REACT_LEVEL, lower_react_level):
            for nbr_atom in atom.get_neighbors():
                if nbr_atom.get_property(ATOM_PROP.K_REACT_LEVEL) is not None:
                    continue
                if with_aromatic and nbr_atom.is_aromatic:
                    cls._sign_react_level_in_aromatic_cycles(nbr_atom, react_level)
                else:
                    nbr_atom.set_property(ATOM_PROP.K_REACT_LEVEL, react_level)

    @classmethod
    def _sign_react_level_in_aromatic_cycles(cls, atom: Atom, react_level: int) -> None:
        for aromatic_cycle in atom.get_aromatic_cycles():
            for atom in aromatic_cycle.get_atoms():
                if atom.get_property(ATOM_PROP.K_REACT_LEVEL) is None:
                    atom.set_property(ATOM_PROP.K_REACT_LEVEL, react_level)

    # ============================
    # -- detect aromatic cycles --

    @classmethod
    def detect_aromatic_cycles(cls, reaction: Reaction) -> None:
        for molecule in reaction.get_molecules():
            AromaticCycleDetector.detect_aromatic_cycles(molecule)

    # ==============================
    # -- check available reaction --
    @classmethod
    def check_available_rxn(cls, reaction: Reaction):
        product_map_num_to_count = cls.check_products_with_full_map(reaction)
        reactant_to_map_nums = cls.check_reactants_with_unique_map_num(reaction)
        # cls.check_one_step_rxn_by_map_num(product_map_num_to_count, reactant_to_map_nums)  # TODO 目前不考虑产物中有多个重复 map num的情况

    # =============================
    # -- check one step reaction --

    @classmethod
    def check_one_step_rxn_by_map_num(cls, product_map_num_to_count: Dict[int, int],
                                      reactant_to_map_nums: Dict[Molecule, List[int]]) -> Dict[Molecule, int]:
        """ 检查是否是一步反应
        基于所有产物中map_num的重复出现次数，
        对每一个反应物，检查其每一个原子的map_num，要求同一个反应物中所有的map_num在产物中的重复出现次数应该相同，否则认为不是单步反应
        TODO 还应该保证产物中相同map_num的原子的react_level相同

        :param product_map_num_to_count:
        :param reactant_to_map_nums:
        :return:
        """
        reactant_to_count = {}
        for reactant, map_nums in reactant_to_map_nums.items():
            count = None
            for map_num in map_nums:
                if count is None:
                    count = product_map_num_to_count[map_num]
                elif count != product_map_num_to_count[map_num]:
                    raise OneStepRxnError('not one step reaction by map num')
            reactant_to_count[reactant] = count
        return reactant_to_count

    @classmethod
    def check_one_step_rxn_by_react_level_and_prop(cls, reaction: Reaction) -> None:
        map_num_to_react_level_and_prop: Dict[int, Tuple[object, int, int]] = {}
        for product in reaction.get_products():
            for atom in product.get_atoms():
                if atom.map_num not in map_num_to_react_level_and_prop.keys():
                    map_num_to_react_level_and_prop[atom.map_num] = (atom.get_property(ATOM_PROP.K_REACT_LEVEL),
                                                                     atom.charge,
                                                                     atom.get_valence())
                else:
                    level_charge_valence = map_num_to_react_level_and_prop[atom.map_num]
                    if level_charge_valence[0] != atom.get_property(ATOM_PROP.K_REACT_LEVEL)\
                            or level_charge_valence[1] != atom.charge\
                            or level_charge_valence[2] != atom.get_valence():
                        raise OneStepRxnError('not one step by react level and charge and valence')

    # =========================================
    # -- check reactants with unique map num --

    @classmethod
    def check_reactants_with_unique_map_num(cls, reaction: Reaction) -> Dict[Molecule, List[int]]:
        passed_map_nums = []
        reactant_to_map_nums = {}
        for reactant in reaction.get_reactants():
            map_nums = []
            for atom in reactant.get_atoms():
                if atom.map_num in passed_map_nums:
                    raise ReactantUniqueMapNumError(f'map num {atom.map_num} is duplicate in reactant')
                elif atom.map_num is not None and atom.map_num != 0:
                    passed_map_nums.append(atom.map_num)
                    map_nums.append(atom.map_num)
            reactant_to_map_nums[reactant] = map_nums
        return reactant_to_map_nums

    # ==================================
    # -- check products with full map --

    @classmethod
    def check_products_with_full_map(cls, reaction: Reaction) -> Dict[int, int]:
        map_num_to_count = {}
        for product in reaction.get_products():
            for atom in product.get_atoms():
                if atom.map_num is None or atom.map_num == 0:
                    raise ProductFullyMapError('product not fully map, some atom has no map num')
                elif not cls._is_reactant_contain_atom_with_map_num(reaction, atom.map_num):
                    raise ProductFullyMapError('product not fully map, some atom not map to reactant')
                if atom.map_num not in map_num_to_count.keys():
                    map_num_to_count[atom.map_num] = 1
                else:
                    raise ReactantDuplicateError('reactant`s equivalent > 1')  # TODO 此处可以通过拆分反应解决
                    map_num_to_count[atom.map_num] += 1
        return map_num_to_count

    @classmethod
    def _is_reactant_contain_atom_with_map_num(cls, reaction: Reaction, map_num: int) -> bool:
        for reactant in reaction.get_reactants():
            for atom in reactant.get_atoms():
                if atom.map_num == map_num:
                    return True
        return False
