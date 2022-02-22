# -*- coding: utf-8 -*-
# @Time     : 2021/4/24 21:28
# @Author   : zbc
# @Email    : zhangbc0315@outlook.com
# @FileName : rxn_center_detector.py
# @Software : PyCharm

from typing import Iterator

from rdkit.Chem import AllChem
from rdkit.Chem.AllChem import ChemicalReaction


class RxnCenterDetector:

    @classmethod
    def _get_atoms_in_mols(cls, mols: [AllChem.Mol]) -> Iterator:
        """ 得到分子列表中的所有原子

        :param mols:
        :return:
        """
        for n, mol in enumerate(mols):
            for atom in mol.GetAtoms():
                yield n, atom

    @classmethod
    def _is_contain_one_atom_by_map_num_in_reactants(cls, reactants: [AllChem.Mol], map_num: int) -> bool:
        """ 反应物列表中是否含有且仅含有一个原子，该原子的map num符合要求的值

        :param reactants:
        :param map_num:
        :return:
        """
        num = 0
        for _, atom in cls._get_atoms_in_mols(reactants):
            if atom.GetAtomMapNum() == map_num:
                num += 1
                if num >= 2:
                    return False
        return num == 1

    @classmethod
    def _check_rxn(cls, rxn: ChemicalReaction) -> None:
        """ 检查反应是否适合抽取反应模板
        1. 反应仅有一个产物
        2. 产物中原子的map num不得相同
        3. 产物中每个原子都有不为0的map num
        4. 产物中原子的原子，在反应物中有且仅有一个原子与之对应

        :param rxn:
        :return:
        """
        if len(rxn.GetProducts()) != 1:
            raise ValueError(f"探测反应中心时，要求有且仅有一个产物，"
                             f"但在{AllChem.ReactionToSmiles(rxn)}中发现了{len(rxn.GetProducts(rxn))}个")
        p_map_nums = set()
        for _, p_atom in cls._get_atoms_in_mols(rxn.GetProducts()):
            p_map_num = p_atom.GetAtomMapNum()
            if p_map_num in p_map_nums:
                raise ValueError(f"探测反应中心时，要求产物中每个原子的map num 不相同，"
                                 f"但在{AllChem.ReactionToSmiles(rxn)}中发现了相同的map num")
            if p_map_num == 0:
                raise ValueError(f"探测反应中心时，要求产物中每个原子都有非0 map num，"
                                 f"但在{AllChem.ReactionToSmiles(rxn)}中发现了0 map num")
            if not cls._is_contain_one_atom_by_map_num_in_reactants(rxn.GetReactants(), p_map_num):
                raise ValueError(f"探测反应中心时，对于产物中的每一个原子，反应物中有且仅有一个原子，与产物中原子的map num相同， "
                                 f"但在{AllChem.ReactionToSmiles(rxn)}中发现错误")
            p_map_nums.add(p_map_num)

    @classmethod
    def _get_atom_in_products_by_map_num(cls, rxn: ChemicalReaction, map_num: int) -> (AllChem.Mol, AllChem.Atom):
        for n, product in enumerate(rxn.GetProducts()):
            for atom in product.GetAtoms():
                if atom.GetAtomMapNum() == map_num:
                    return n, atom
        return n, None

    @classmethod
    def _dict_append(cls, dic: {}, mol: AllChem.Mol, val: int):
        if mol not in dic.keys():
            dic[mol] = []
        dic[mol].append(val)

    @classmethod
    def _is_atom_changed(cls, r_atom, p_atom) -> bool:
        """ 判断反应物原子到产物原子是否发生了变化
        1. 如果原子的临近原子数目发生变化，则 True
        2. 对于反应物原子的任一临近原子，无法在产物原子的临近原子中找到相同 map num 的原子，则 True

        :param r_atom:
        :param p_atom:
        :return:
        """
        if p_atom is None:
            return True
        r_nbr_atoms = r_atom.GetNeighbors()
        p_nbr_atoms = p_atom.GetNeighbors()
        r_nbr_atoms_map_nums = [ar.GetAtomMapNum() for ar in r_nbr_atoms]
        p_nbr_atoms_map_nums = [ap.GetAtomMapNum() for ap in p_nbr_atoms]
        r_nbr_atoms_idxes = [ar.GetIdx() for ar in r_nbr_atoms]
        p_nbr_atoms_idxes = [ap.GetIdx() for ap in p_nbr_atoms]
        if len(r_nbr_atoms_map_nums) != len(p_nbr_atoms_map_nums):
            return True
        for r_nbr_atom_map_num, r_nbr_atom_idx in zip(r_nbr_atoms_map_nums, r_nbr_atoms_idxes):
            if r_nbr_atom_map_num not in p_nbr_atoms_map_nums:
                return True
            p_nbr_atom_idx = p_nbr_atoms_idxes[p_nbr_atoms_map_nums.index(r_nbr_atom_map_num)]
            r_bond = r_atom.GetOwningMol().GetBondBetweenAtoms(r_atom.GetIdx(), r_nbr_atom_idx)
            p_bond = p_atom.GetOwningMol().GetBondBetweenAtoms(p_atom.GetIdx(), p_nbr_atom_idx)
            if r_bond.GetBondType() != p_bond.GetBondType():
                return True
        return False

    @classmethod
    def detect_rxn_center(cls, rxn: ChemicalReaction) -> ({AllChem.Mol}, {AllChem.Mol}):
        cls._check_rxn(rxn)
        reactants_center = {}
        product_center = {}
        for reactant_i, r_atom in cls._get_atoms_in_mols(rxn.GetReactants()):
            r_map_num = r_atom.GetAtomMapNum()
            r_idx = r_atom.GetIdx()
            if r_map_num == 0:
                cls._dict_append(reactants_center, reactant_i, r_idx)
                continue
            product_i, p_atom = cls._get_atom_in_products_by_map_num(rxn, r_map_num)
            if cls._is_atom_changed(r_atom, p_atom):
                cls._dict_append(reactants_center, reactant_i, r_idx)
                if p_atom is not None:
                    cls._dict_append(product_center, product_i, p_atom.GetIdx())
        return reactants_center, product_center


if __name__ == '__main__':
    pass
