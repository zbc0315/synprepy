# -*- coding: utf-8 -*-
# @Time     : 2020/4/20 17:16
# @Author   : zbc
# @Email    : zhangbc0315@outlook.com
# @FileName : smiles_parser.py
# @Software : PyCharm

from typing import List

from crkitpy.smilesparse.smiles_lex import SmilesLex
from crkitpy.smilesparse.smiles_yacc import SmilesYacc
from crkitpy.molgraph.reaction import Reaction
from crkitpy.molgraph.molecule import Molecule
from crkitpy.molgraph.atom import Atom, Hybrid
from crkitpy.smilesparse.rdkit_smiles_parser import RdkitSmilesParser
from crkitpy.chemdata.atomic_data import AtomicData
from crkitpy.molhandler.cycle_detector import CycleDetector
from crkitpy.const import BOND_TYPE


class SmilesParser:

    # ========================================
    # -- parse to molecule default by rdkit --

    @classmethod
    def parse_to_mol(cls, smiles: str) -> Molecule:
        return RdkitSmilesParser.parse_to_mol(smiles)
        # smiles_lex = SmilesLex()
        # smiles_yacc = SmilesYacc()
        # molecule = smiles_yacc.parse(smiles_lex.input(smiles))
        # return molecule

    # ==============================
    # -- parse to molecule by lex --

    @classmethod
    def parse_to_mol_by_lex(cls, smiles: str) -> Molecule:
        smiles_lex = SmilesLex()
        smiles_yacc = SmilesYacc()
        molecule = smiles_yacc.parse(smiles_lex.input(smiles))
        cls._set_valence_and_hybrid(molecule)
        cls._get_aromatic_hybrid_by_hybrid(molecule)
        return molecule

    # =================================
    # -- get aromatic info by hybrid --

    @classmethod
    def _get_aromatic_hybrid_by_hybrid(cls, molecule: Molecule) -> None:
        # TODO 此处认为环中所有原子皆为SP2杂化或可以为SP2杂化时，设定环为芳香环，
        #      不排除有些更复杂的杂化类型也可以出现在芳香环中
        #      不排除环中全是SP2杂化，但是环依然不具有芳香性
        for cycle in CycleDetector.get_cycles(molecule):
            if all([atom.get_hybrid() == Hybrid.SP2 or Hybrid.SP2 in atom.get_alter_hybrids() for atom in cycle]):
                for atom in cycle:
                    atom.set_hybrid(Hybrid.SP2)
                    atom.is_aromatic = True
                    for bond in atom.get_bonds():
                        if any([a not in cycle for a in bond.get_atoms()]):
                            continue
                        bond.set_bond_type(BOND_TYPE.AROMATIC)

    # ======================================
    # -- correct hybrid by aromatic cycle --

    # @classmethod
    # def _correct_hybrid_by_aromatic_cycle(cls, molecule: Molecule) -> None:
    #     """ 根据芳香环信息更正原子的杂化情况
    #     例如，呋喃中的氧从SP3更正为SP2
    #
    #     """
    #     for cycle in CycleDetector.get_cycles(molecule):
    #         if all([atom.get_hybrid() == Hybrid.SP2 or Hybrid.SP2 in atom.get_alter_hybrids() for atom in cycle]):
    #             for atom in cycle:
    #                 atom.set_hybrid(Hybrid.SP2)
    #                 atom.is_aromatic = True

    # ========================
    # -- valence and hybrid --

    @classmethod
    def _set_valence_and_hybrid(cls, molecule: Molecule) -> None:
        """ 此处应根据每个原子与周围非氢原子的成键类型设置本原子的explicit valence
        由于芳香键对不同类型原子在不同情况下所代表的valence不同：
            例如对于连接两个芳香键的N原子，如果其还连接一个H原子，则其explicit valence为2
                                      如果其不连接H原子，则其explicit valence为3
        TODO 设置valence字典，根据字典设置每个原子的valence
        :param molecule:
        :return:
        """
        for atom in molecule.get_atoms():
            if atom.is_r_group():
                continue
            atomic_data = AtomicData.get_element_data(atom.get_symbol())
            cls._set_atom_explicit_valence(atom, atomic_data.valences)
            cls._set_atom_hybrid(atom, atomic_data.outshell_elec_num)

    @classmethod
    def _set_atom_explicit_valence(cls, atom: Atom, valences: List[int]) -> None:
        bond_weight = atom.get_bonds_weight_without_h()
        atom.set_explicit_valence(bond_weight)  # TODO 此处未处理隐式氢
        # TODO 此处未处理电荷引起的价态变化
        """
        for valence in valences:
            if valence >= bond_weight + atom.get_hs_num():
                atom.set_implicit_valence(valence - bond_weight)  
                break
        """

    @classmethod
    def _set_atom_hybrid(cls, atom: Atom, oe_num: int) -> None:
        """
        :param atom: 原子
        :param oe_num: 原子最外层电子数
        :return:
        """
        # 未参与成键的最外层电子数
        roe_num = oe_num - atom.get_valence()
        if roe_num < 0:
            roe_num = 0  # TODO 例如P最外层有5个电子，但却可以有7价
        roe_pair_num = int(roe_num/2)
        hybrid = Hybrid(roe_pair_num + atom.get_hs_num() + len(atom.get_bonds()))
        atom.set_hybrid(hybrid)
        if roe_pair_num > 0 and atom.get_hybrid() == Hybrid.SP3:
            atom.add_alter_hybrid(Hybrid(atom.get_hybrid().value - 1))

    # =======================
    # -- parse to reaction --

    @classmethod
    def parse_to_rxn(cls, smiles: str) -> Reaction:
        return RdkitSmilesParser.parse_to_rxn(smiles)


if __name__ == '__main__':
    smi = '[CH3:25][CH:21]([NH:20][P:19](=[O:26])([O:18][CH2:17][CH:14]1[O:13][C:8]([C:9]#[N:10])([CH:11]([OH:12])[CH:15]1[OH:16])[c:7]1[cH:34][cH:35][c:36]2[c:2]([NH2:1])[nH:3][cH:4][nH:5][n:6]12)[O:27][c:28]1[cH:29][cH:30][cH:31][cH:32][cH:33]1)[C:22]([OH:24])=[O:23]'
    mol = SmilesParser.parse_to_mol_by_lex(smi)
    mol.draw()
