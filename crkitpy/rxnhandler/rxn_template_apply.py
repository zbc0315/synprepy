# -*- coding: utf-8 -*-
# @Time     : 2020/6/4 15:51
# @Author   : zbc
# @Email    : zhangbc0315@outlook.com
# @FileName : rxn_template_apply.py
# @Software : PyCharm

from typing import List, Iterator, Dict, Tuple, Union

from memory_profiler import profile
import networkx as nx
from networkx.algorithms import isomorphism

from crkitpy.molhandler.cycle_detector import CycleDetector
from crkitpy.molgraph.reaction import Reaction, Molecule
from crkitpy.molgraph.atom import Atom
from crkitpy.const import BOND_TYPE, ATOM_PROP
from crkitpy.const.bond_type import BOND_WEIGHTS_TO_TYPE
from crkitpy.errors.rxn_error import TemplateGroupDuplicateError
from crkitpy.errors.bond_error import BondError
from crkitpy.rxnhandler.reset_hs_num_after_apply import ResetHsNumAfterApply
from crkitpy.smilesparse import SmilesParser
from crkitpy.smatesparse.smates_parser import SmatesParser


class RxnTemplateApply:

    """
    0: None
    1: [NH3:1].[OH:2][CH3:3]>>[NH2:1][CH3:3]
    """
    _temp_type = 0

    # ===========================
    # -- reactants to products --

    @classmethod
    def reactants_to_products(cls, smates: str, reactants: List[Molecule], strict: bool=False) -> Iterator[List[Molecule]]:
        cls._temp_type = 0
        if smates == '[NH3:1].[OH:2][CH3:3]>>[NH2:1][CH3:3]':
            cls._temp_type = 1
        num_reactants = len(reactants)
        union_reactant = Molecule.union_mols(reactants)
        reverse_smates = cls._reverse_smates(smates)
        # rxn_template = cls._reverse_rxn_template(rxn_template)
        for products_and_num_joined in cls.product_to_reactants(reverse_smates, union_reactant, True, need_num_joined_mol=True):
            if not strict or num_reactants == products_and_num_joined[-1]:
                yield products_and_num_joined[0]

    # ==========================
    # -- reverse rxn template --

    @classmethod
    def _reverse_smates(cls, smates: str) -> str:
        smates_r, smates_p = smates.split('>>')
        return f'{smates_p}>>{smates_r}'

    @classmethod
    def _reverse_rxn_template(cls, rxn_template: Reaction) -> Reaction:
        res = Reaction()
        res.add_product(rxn_template.get_reactants()[0])
        res.add_reactant(rxn_template.get_products()[0])
        return res

    # ==========================
    # -- product to reactants --

    @classmethod
    def product_to_reactants(cls, smates: str, product: Molecule,
                             reverse: bool = False, need_num_joined_mol: bool = False) -> Iterator[Union[List[Molecule], Tuple[List[Molecule], int]]]:
        """ 使用rxn template从产物预测得到可能的反应物

        :param reverse:
        :param smates:
        :param product:
        :param need_num_joined_mol:
        :return: 反应物集合的生成器，因为在将模板对产物进行子图同构映射时，可能会产生多种映射方案，所以有多组可能的反应物集合
        """
        cls._temp_type = 0
        if smates in ['[NH3:1].[OH:2][CH3:3]>>[NH2:1][CH3:3]', '[NH2:1][CH3:3]>>[NH3:1].[OH:2][CH3:3]']:
            cls._temp_type = 1
        res = []
        num_joined_mol = []
        rxn_template = SmatesParser.template_from_smates(smates)
        for mapping_t_p_atoms in cls._match_product_template(rxn_template.get_products()[0], product):
            # cls._reset_hs_num(mapping_t_p_atoms, rxn_template, reverse)
            reactants = cls._predict_reactants(mapping_t_p_atoms, rxn_template.get_reactants()[0].copy()[0], product)
            try:
                ResetHsNumAfterApply.reset_hs_num(reactants)
            except Exception as e:
                continue
            # res.append(reactants)
            if need_num_joined_mol:
                res.append((reactants, cls._count_num_joined_mol(mapping_t_p_atoms)))
            else:
                res.append(reactants)
            # yield reactants
        if need_num_joined_mol:
            res = sorted(res, key=lambda x: len(x[0]), reverse=True)
        else:
            res = sorted(res, key=lambda x: len(x), reverse=True)
        for rs in res:
            yield rs

    @classmethod
    def _count_num_joined_mol(cls, mapping_t_p_atoms):
        mids = set()
        for p_atom in mapping_t_p_atoms.values():
            mids.add(p_atom.get_property(ATOM_PROP.K_MOL_ID))
        # mids.remove(None)
        return len(mids)

    # ===================================

    @classmethod
    def _reset_hs_num(cls, mapping_t_p_atoms, rxn_template, reverse: bool):
        p = rxn_template.get_products()[0] if not reverse else rxn_template.get_reactants()[0]
        r = rxn_template.get_reactants()[0] if not reverse else rxn_template.get_products()[0]
        for t_atom in mapping_t_p_atoms.keys():
            p_atom = mapping_t_p_atoms[t_atom]
            tr_atom = cls._get_atom_by_map_num_in_mols(r, t_atom.map_num)
            hs_num = tr_atom.get_hs_num()
            tr_atom.explicit_hs_num = 0
            new_hs_num = hs_num + (p_atom.get_hs_num() - t_atom.get_hs_num())
            tr_atom.implicit_hs_num = new_hs_num if new_hs_num >= 0 else 0

    @classmethod
    def _get_atom_by_map_num_in_mols(cls, mol, map_num):
        atom = mol.get_one_atom_by_map_num(map_num)
        return atom

    @classmethod
    def _predict_reactants(cls, mapping_t_p_atoms: Dict[Atom, Atom], reactants_template: Molecule, product: Molecule):
        duplicate_r_group_atoms = cls._mark_and_extract_duplicate_group(reactants_template, mapping_t_p_atoms)
        save_p_atoms = []
        for p_atom in product.get_atoms():
            if p_atom not in mapping_t_p_atoms.values():
                save_p_atoms.append(p_atom)
        save_product = product.copy_sub_mol(save_p_atoms)

        if len(duplicate_r_group_atoms) > 0:
            dup_all_groups = []
            dup_all_link_this_other_atoms = []
            try:
                for dup_group, dup_link_this_other_atoms in cls._copy_duplicate_groups(save_product):
                    dup_all_groups.append(dup_group)
                    dup_all_link_this_other_atoms.extend(dup_link_this_other_atoms)
            except TemplateGroupDuplicateError as e:
                print(e)
                return []
            union_dup_group = Molecule.union_mols(dup_all_groups)
            reactants_template = reactants_template.add_molecule(union_dup_group, dup_all_link_this_other_atoms)

        link_this_other_atoms = []
        for tp_atom, p_atom in mapping_t_p_atoms.items():
            # if tp_atom in duplicate_r_group_atoms:
            #     continue
            # if tp_atom.map_num == 13:
            #     print(13)
            tr_atom = reactants_template.get_one_atom_by_map_num(tp_atom.map_num,
                                                                 lambda x: not x.get_property(
                                                                     ATOM_PROP.K_DO_NOT_LINK_OTHERS))
            if tr_atom is None:
                # TODO 舍弃的部分也应补全
                continue
            if tr_atom.is_r_group():
                tr_atom.copy_atom(p_atom, without_map_num=True)
            for p_nbr_atom in p_atom.get_neighbors():
                if p_nbr_atom in mapping_t_p_atoms.values():
                    continue
                bond = p_atom.get_bond(p_nbr_atom)
                # tr_atom.implicit_hs_num = tr_atom.implicit_hs_num - cls._bond_type_to_num(bond.get_bond_type())
                # if tr_atom.implicit_hs_num < 0:
                #     tr_atom.implicit_hs_num = 0
                tr_atom.set_property(ATOM_PROP.K_REACT_CENTER, True)
                link_this_other_atoms.append((save_product.get_atom_by_parent(p_nbr_atom),
                                              bond.get_bond_type(),
                                              tr_atom))
        res = save_product.add_molecule(reactants_template, link_this_other_atoms)
        return cls._get_right_mols(res.get_detached_sub_mols())

    @classmethod
    def _get_right_mols(cls, mols):
        res = []
        for mol in mols:
            if cls._is_mapped_atom_in_mol(mol):
                res.append(mol)
        return res

    @classmethod
    def _is_mapped_atom_in_mol(cls, mol):
        for atom in mol.get_atoms():
            if atom.map_num is not None and atom.map_num != 0:
                return True
        return False

    @classmethod
    def _bond_type_to_num(cls, bond_type):
        d = {BOND_TYPE.SINGLE: 1,
             BOND_TYPE.DOUBLE: 2,
             BOND_TYPE.TRIPLE: 3,
             BOND_TYPE.AROMATIC: 1}
        return d[bond_type]

    # ==========================
    # -- copy duplicate group --

    @classmethod
    def _copy_duplicate_groups(cls, save_product: Molecule) -> Iterator[Tuple[Molecule, List[Tuple[Atom, str, Atom]]]]:
        for group in save_product.get_detached_sub_mols():
            for new_group, link_this_other_atoms in cls._copy_duplicate_group(group):
                yield new_group, link_this_other_atoms

    @classmethod
    def _copy_duplicate_group(cls, group: Molecule):
        group_duplicate_info = None
        tobe_connected_atoms = []
        for atom in group.get_atoms():
            atom_duplicate_info = atom.get_property(ATOM_PROP.K_PART_DUPLICATE_TIMES)
            if atom_duplicate_info is None:
                continue
            if group_duplicate_info is None:
                group_duplicate_info = atom_duplicate_info
            elif group_duplicate_info != atom_duplicate_info:
                raise TemplateGroupDuplicateError("duplicated group is not separated")
            tobe_connected_atoms.append(atom)
        if group_duplicate_info is not None:
            times = group_duplicate_info[-1]
            p_atom = group_duplicate_info[0]
            tr_atoms = group_duplicate_info[1][:-1]
            link_this_other_atoms = []
            for t in range(times - 1):
                tr_atom = tr_atoms[t]
                new_group, _, __ = group.copy(save_parent=True, save_atom_copy_hash=False)
                for tobe_connected_atom in tobe_connected_atoms:
                    new_atom = new_group.get_atom_by_parent(tobe_connected_atom)
                    p_nbr_atom = tobe_connected_atom.get_parent().get_parent()
                    bond = p_atom.get_bond(p_nbr_atom)
                    link_this_other_atoms.append((new_atom, bond.get_bond_type(), tr_atom))
                tr_atom.set_property(ATOM_PROP.K_DO_NOT_LINK_OTHERS, True)
                yield new_group, link_this_other_atoms

    # ==========================
    # -- mark duplicate group --

    @classmethod
    def _mark_and_extract_duplicate_group(cls, reactants_template: Molecule, mapping_t_p_atoms: Dict[Atom, Atom]) -> \
            List[Atom]:
        """
        :param reactants_template:
        :param mapping_t_p_atoms:
        :return:
        """
        res = []
        for tp_atom, p_atom in mapping_t_p_atoms.items():
            tr_atoms = list(reactants_template.get_atoms_by_map_num(tp_atom.map_num))
            if len(tr_atoms) > 1:
                for nbr_atom in p_atom.get_neighbors():
                    nbr_atom.set_property(ATOM_PROP.K_PART_DUPLICATE_TIMES, (p_atom, tr_atoms, len(tr_atoms)))
                res.append(tp_atom)
                if any([tr_atom.is_r_group() for tr_atom in tr_atoms]):
                    [tr_atom.copy_atom(p_atom, without_map_num=True) for tr_atom in tr_atoms]
        return res

    # ============================
    # -- match product template --

    @classmethod
    def _match_product_template(cls, product_template: Molecule, product: Molecule) -> Iterator[Dict[Atom, Atom]]:
        gm = isomorphism.GraphMatcher(product.get_mol_graph_visible(),
                                      product_template.get_mol_graph_visible(),
                                      node_match=cls._node_match_func,
                                      edge_match=cls._edge_match_func)
        res = gm.subgraph_isomorphisms_iter()
        for mapping in res:
            # product_template.draw()
            # cls._show_match_result(mapping, product)
            yield dict(zip(mapping.values(), mapping.keys()))

    @classmethod
    def _show_match_result(cls, match_p_t_atoms: {Atom: Atom}, product: Molecule):
        for p_atom, t_atom in match_p_t_atoms.items():
            p_atom.set_property(ATOM_PROP.K_NEED_SAVE, True)

    @classmethod
    def _node_match_func(cls, node, t_node) -> bool:
        if cls._temp_type == 1 and node['symbol'] == 'O' and node['hs_num'] == 0:
            return False
        if t_node['symbol'] == '*':
            return True
        elif t_node['symbol'] == 'ARO' and node['is_aromatic']:
            return True
        # print(f"{node['symbol']} - {node['no_h_bond_valence']}, {t_node['symbol']} - {t_node['require_valence']} - {t_node['max_valence']}")
        return node['symbol'] == t_node['symbol'] and \
               node['charge'] == t_node['charge'] and \
               node['is_aromatic'] == t_node['is_aromatic'] and \
               node['no_h_bond_valence'] + t_node['require_valence'] <= t_node['max_valence']

    @staticmethod
    def _edge_match_func(edge, t_edge) -> bool:
        return edge['bond_type'] == t_edge['bond_type']
        # return edge['bond_type'] == t_edge['bond_type'] \
        #        or (edge['bond_type'] == BOND_TYPE.AROMATIC and t_edge['bond_type'] in [BOND_TYPE.DOUBLE,
        #                                                                                BOND_TYPE.SINGLE]) \
        #        or (edge['bond_type'] in [BOND_TYPE.DOUBLE, BOND_TYPE.SINGLE] and t_edge[
        #     'bond_type'] == BOND_TYPE.AROMATIC)

    # ===============================
    # --            drop           --
    # ===============================

    '''
    # ==================================
    # -- count reactant duplicate num --

    @classmethod
    def _count_reactant_duplicate_num(cls, tp_map_num_to_count: Dict[int, int],
                                      tr_to_map_nums: Dict[Molecule, List[int]]) -> Dict[Molecule, int]:
        """ 获得反应模板中，每个反应物的重复次数（当量）

        :param tp_map_num_to_count: 模板产物中map num 和 map num出现次数的键值映射
        :param tr_to_map_nums: 模板反应物和其map nums的键值映射
        :return: 模板中反应物 - 该反应物的重复次数（当量）
        """
        return RxnTemplateExtractor.check_one_step_rxn_by_map_num(tp_map_num_to_count, tr_to_map_nums)

    # =======================================
    # -- count map num in template product --

    @classmethod
    def _count_map_num_in_template_product(cls, rxn_template_selector: Reaction) -> Dict[int, int]:
        """ 统计模板产物中map num的重复次数
        会检查是否每个模板产物中的原子都有map num，如果否，则报错
        会检查是否每个模板产物中的原子都有和它map num相同的模板反应物中的原子，如果否，则报错
        TODO 此处处理忽略了多产物反应模板可能导致的后续错误，如果尝试使用多产物反应模板预测一个产物的可能反应物，只有一个模板产物会与产物
             出现子图同构映射，但是通过本方法得到的是所有模板产物中的map num的count，在后续操作中会出错，因此，当前仅可以使用单产物反应
             模板。
        :param rxn_template_selector:
        :return: 模板产物中的原子的map num - 该map num在模板产物中的重复次数
        """
        return RxnTemplateExtractor.check_products_with_full_map(rxn_template_selector)
    '''


if __name__ == '__main__':
    # pass
    from crkitpy.smatesparse.smates_parser import SmatesParser
    from crkitpy.inchiparse.inchi_parser import InchiParser
    from crkitpy.smilesparse.smiles_writer import SmilesWriter
    from crkitpy.inchiparse.inchi_writer import InchiWriter
    sma = "[Cl:2][CH3:3].[NH3:1]>>[NH2:1][CH3:3]"
    # sma = "[CH4:1].[CH4:2]>>[CH3:1][CH3:2]"
    # rs_inchis = ["InChI=1S/C9H9Cl/c10-9-5-7-3-1-2-4-8(7)6-9/h1-4,9H,5-6H2", "InChI=1S/C4H9N/c5-4-2-1-3-4/h4H,1-3,5H2"]
    rs_smiles = ['O=C(N1CCNCC1)c2ccccc2', 'COc1c2c(C(C(O)=O)=O)cn(COP(O)(O)=O)c2c(n3nc(C)nc3)nc1']
    sma = '[NH3:1].[OH:2][CH3:3]>>[NH2:1][CH3:3]'

    # sma = "[N:1]#[CH:3].[NH2:2][NH2:4]>>[N:1](=[CH2:3])[NH2:2]"
    # rs_inchi = ["InChI=1S/C3H2ClN/c1-3(4)2-5/h1H2", "InChI=1S/C7H5Cl2F3N2/c8-4-1-3(7(10,11)12)2-5(9)6(4)14-13/h1-2,14H,13H2"]

    # reactants = list(map(InchiParser.mol_from_inchi, rs_inchis))
    temp = SmatesParser.template_from_smates(sma)
    rs_smiles = ['CC(OCC)C', 'NC1CCCCC1']
    sma = '[CH3:2][CH2:5][O:4][CH3:3].[NH3:1]>>[NH2:1][CH3:3]'

    reactants = list(map(SmilesParser.parse_to_mol, rs_smiles))
    for products in RxnTemplateApply.reactants_to_products(smates=sma, reactants=reactants, strict=True):
        # print(products[-1])
        # if len(products) > 1:
        for p in products:
            print(InchiWriter.mol_to_inchi(p))
        # print(InchiWriter.mol_to_inchi(products[0]))
        print("="*20)
