#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Time    : 2022/1/5 10:43
# @Author  : zhangbc0315@outlook.com
# @File    : template_extractor.py
# @Software: PyCharm

from enum import Enum

from rdkit.Chem.rdchem import Atom, Bond, Mol
from rdkit.Chem import AllChem

from chem_utils.chem_handler import AtomHandler
from chem_utils.chem_handler import MolHandler
from chem_utils.chem_handler import RxnHandler


class RxnTemplateType(Enum):

    CENTRALIZED = 0,
    EXTENDED = 1


class TemplateExtractor:

    # region ===== assign neighbor atoms =====

    @classmethod
    def _assign_atom_nbr_as_nbr(cls, atom: Atom):
        for nbr in atom.GetNeighbors():
            if AtomHandler.is_atom_reacted(nbr) or AtomHandler.is_atom_dropped(nbr):
                continue
            else:
                AtomHandler.set_atom_is_nbr(nbr)

    @classmethod
    def _assign_rxn_center_nbr_atoms_to_mols(cls, mols: [Mol]):
        for mol in mols:
            for atom in mol.GetAtoms():
                if AtomHandler.is_atom_reacted(atom):
                    cls._assign_atom_nbr_as_nbr(atom)

    @classmethod
    def _assign_rxn_center_nbr_atoms(cls, rxn: AllChem.ChemicalReaction):
        cls._assign_rxn_center_nbr_atoms_to_mols(rxn.GetReactants())
        cls._assign_rxn_center_nbr_atoms_to_mols(rxn.GetProducts())

    # endregion

    # region ===== assign special atoms =====

    @classmethod
    def _is_atom_dropped(cls, map_num: int, rxn: AllChem.ChemicalReaction) -> bool:
        return map_num is None or map_num == 0 \
               or not MolHandler.contain_atom_with_map_num_in_mols(rxn.GetProducts(), map_num)

    @classmethod
    def _has_same_num_neighbor(cls, r_atom: Atom, p_atom: Atom) -> bool:
        return len(r_atom.GetNeighbors()) == len(p_atom.GetNeighbors())

    @classmethod
    def _get_nbr_with_map_num(cls, atom: Atom, map_num: int) -> Atom:
        for nbr_atom in atom.GetNeighbors():
            if nbr_atom.GetAtomMapNum() == map_num:
                return nbr_atom
        return None

    @classmethod
    def _is_same_atoms(cls, atom1: Atom, atom2: Atom) -> bool:
        return atom1.GetAtomicNum() == atom2.GetAtomicNum() \
               and atom1.GetNumRadicalElectrons() == atom2.GetNumRadicalElectrons() \
               and len(atom1.GetNeighbors()) == len(atom2.GetNeighbors())

    @classmethod
    def _is_same_bond(cls, bond1: Bond, bond2: Bond):
        return bond1.GetBondType() == bond2.GetBondType()

    @classmethod
    def _has_same_neighbor_atom_with_same_bond(cls, r_atom: Atom, p_atom: Atom) -> bool:
        for r_nbr_atom in r_atom.GetNeighbors():
            p_nbr_atom = cls._get_nbr_with_map_num(p_atom, r_nbr_atom.GetAtomMapNum())
            if p_nbr_atom is None:
                return False
            r_bond = r_atom.GetOwningMol().GetBondBetweenAtoms(r_atom.GetIdx(), r_nbr_atom.GetIdx())
            p_bond = p_atom.GetOwningMol().GetBondBetweenAtoms(p_atom.GetIdx(), p_nbr_atom.GetIdx())
            if not cls._is_same_bond(r_bond, p_bond):
                return False
        return True

    @classmethod
    def _is_atom_reacted(cls, r_atom: Atom, p_atom: Atom) -> bool:
        if not cls._is_same_atoms(r_atom, p_atom):
            return True
        if not cls._has_same_neighbor_atom_with_same_bond(r_atom, p_atom):
            return True
        return False

    @classmethod
    def _assign_dropped_and_reacted_atoms_in_rxn(cls, rxn: AllChem.ChemicalReaction):
        for reactant in rxn.GetReactants():
            for atom in reactant.GetAtoms():
                map_num = atom.GetAtomMapNum()
                if cls._is_atom_dropped(map_num, rxn):
                    AtomHandler.set_atom_is_dropped(atom)
                    atom.SetAtomMapNum(0)
                else:
                    p_atom = list(MolHandler.query_atoms_with_map_num_in_mols(rxn.GetProducts(), map_num))[0]
                    if cls._is_atom_reacted(atom, p_atom):
                        AtomHandler.set_atom_is_reacted(atom)
                        AtomHandler.set_atom_is_reacted(p_atom)

    # endregion

    # region ===== check =====

    @classmethod
    def _check_map_num_unique_in_mols(cls, mols: [Mol]):
        checked_map_num = set()
        for mol in mols:
            for atom in mol.GetAtoms():
                map_num = atom.GetAtomMapNum()
                if map_num is None or map_num == 0:
                    continue
                if map_num in checked_map_num:
                    raise ValueError(f"There are more than one atoms have same the map num: {map_num} in reactants or products")
                checked_map_num.add(map_num)

    @classmethod
    def _check_map_num_unique(cls, reaction: AllChem.ChemicalReaction):
        cls._check_map_num_unique_in_mols(reaction.GetReactants())
        cls._check_map_num_unique_in_mols(reaction.GetProducts())

    @classmethod
    def _check_map_num_enough(cls, reaction: AllChem.ChemicalReaction):
        for product in reaction.GetProducts():
            for atom in product.GetAtoms():
                map_num = atom.GetAtomMapNum()
                if not MolHandler.contain_atom_with_map_num_in_mols(reaction.GetReactants(), map_num):
                    raise ValueError(f"Can not find atom in reactants contain map num({map_num}) same with atom in product")

    @classmethod
    def _check_rxn_contain_reactant_and_product(cls, rxn: AllChem.ChemicalReaction):
        if rxn.GetNumProductTemplates() == 0 or rxn.GetNumReactantTemplates() == 0:
            raise ValueError(f"The reaction or template don't contain reactants or products")

    @classmethod
    def _check_rxn_contain_one_product(cls, rxn: AllChem.ChemicalReaction):
        if rxn.GetNumProductTemplates() != 1:
            raise ValueError(f"Only need one product, but get: {rxn.GetNumProductTemplates()}")

    # endregion

    # region ===== get reaction center =====

    @classmethod
    def _should_atom_saved(cls, atom: Atom) -> bool:
        return AtomHandler.is_atom_dropped(atom) or AtomHandler.is_atom_reacted(atom) or AtomHandler.is_atom_nbr(atom)

    @classmethod
    def _get_mol_template(cls, mol: Mol) -> Mol:
        atoms_should_saved_idxes = []
        # print(AllChem.MolToSmiles(mol))
        for atom in mol.GetAtoms():
            if cls._should_atom_saved(atom):
                atoms_should_saved_idxes.append(atom.GetIdx())
        if len(atoms_should_saved_idxes) > 0:
            return MolHandler.frag_idxes_to_mol(mol, atoms_should_saved_idxes)
        else:
            return None

    @classmethod
    def _get_rxn_template(cls, rxn: AllChem.ChemicalReaction) -> AllChem.ChemicalReaction:
        rxn_template = AllChem.ChemicalReaction()
        reactants = []
        for reactant in rxn.GetReactants():
            reactant_template = cls._get_mol_template(reactant)
            if reactant_template is not None:
                reactants.append(reactant_template)
        if len(reactants) > 0:
            rxn_template.AddReactantTemplate(MolHandler.merge_mols(reactants))
        products = []
        for product in rxn.GetProducts():
            product_template = cls._get_mol_template(product)
            if product_template is not None:
                products.append(product_template)
        if len(products) > 0:
            rxn_template.AddProductTemplate(MolHandler.merge_mols(products))
        return rxn_template

    # endregion =====

    @classmethod
    def canon_rxn_template(cls, rxn_template: AllChem.ChemicalReaction):
        map_num_old_to_new = MolHandler.canonicalize_map_num(rxn_template.GetProducts()[0])
        for atom in rxn_template.GetReactants()[0].GetAtoms():
            map_num = atom.GetAtomMapNum()
            if map_num != 0:
                atom.SetAtomMapNum(map_num_old_to_new[map_num])
        AllChem.CanonicalizeMol(rxn_template.GetReactants()[0])

    @classmethod
    def rxn_smiles_to_rxn_temp_smarts(cls, rxn_smiles: str, rxn_template_type: RxnTemplateType, need_clean: bool) -> str:
        rxn = RxnHandler.smarts_to_rxn(rxn_smiles, True)
        if need_clean:
            rxn = RxnHandler.remove_unmapped_mols_in_rxn(rxn)
            rxn = RxnHandler.remove_products_same_with_reactants(rxn)
        cls._check_rxn_contain_reactant_and_product(rxn)
        cls._check_rxn_contain_one_product(rxn)
        cls._check_map_num_unique(rxn)
        cls._check_map_num_enough(rxn)
        cls._assign_dropped_and_reacted_atoms_in_rxn(rxn)
        if rxn_template_type == RxnTemplateType.EXTENDED:
            cls._assign_rxn_center_nbr_atoms(rxn)
        rxn.Initialize()
        rxn_template = cls._get_rxn_template(rxn)
        cls._check_rxn_contain_reactant_and_product(rxn_template)
        cls.canon_rxn_template(rxn_template)
        return AllChem.ReactionToSmiles(rxn_template)


if __name__ == "__main__":

    # smi = "CS[c:31]1[n:30]-[c:29]2[c:34]([cH:33][n:32]1)[C:35](=[O:36])[N:26]([c:20]1[c:19]([Cl:18])[cH:24][cH:23][cH:22][c:21]1[Cl:25])[CH2:27][N:28]2[CH2:39][CH3:40].[NH2:1][c:2]1[cH:3][cH:4][c:5]([N:10]2[CH2:11][C@H:12]([CH3:17])[NH:13][C@H:14]([CH3:16])[CH2:15]2)[c:6]([CH2:8][OH:9])[cH:7]1>>[NH:1]([c:2]1[cH:3][cH:4][c:5]([N:10]2[CH2:11][C@H:12]([CH3:17])[NH:13][C@H:14]([CH3:16])[CH2:15]2)[c:6]([CH2:8][OH:9])[cH:7]1)[c:31]1[n:30]-[c:29]2[c:34]([cH:33][n:32]1)[C:35](=[O:36])[N:26]([c:20]1[c:19]([Cl:18])[cH:24][cH:23][cH:22][c:21]1[Cl:25])[CH2:27][N:28]2[CH2:39][CH3:40]"
    smi = "O=S(Cl)[Cl:3].O[CH3:21].[NH2:5][CH:6]([C:7](=[O:8])[OH:9])[C:10]1=[CH:11][CH:12]=[C:13]([O:16][C:17]([F:18])([F:19])[F:20])[CH:14]=[CH:15]1>>[ClH:3].[NH2:5][CH:6]([C:7](=[O:8])[O:9][CH3:21])[C:10]1=[CH:15][CH:14]=[C:13]([O:16][C:17]([F:18])([F:19])[F:20])[CH:12]=[CH:11]1"
    # smi = "[CH2:1]([C:2]1=[CH:3][CH:4]=[CH:5][CH:6]=[CH:7]1)[CH:8]1[O:9][CH2:10]1>>[CH2:1]([C:2]1=[CH:3][CH:4]=[CH:5][CH:6]=[CH:7]1)[C@@H:8]1[O:9][CH2:10]1"
    # smi = "[ClH:1].[N:2]1([CH:11]([CH:12]([CH2:13][NH:14][CH3:15])[OH:16])[C:17]2=[CH:18][CH:19]=[C:20]([O:23][C:24]([F:25])([F:26])[F:27])[CH:21]=[CH:22]2)[CH:3]=[CH:4][C:5]2=[CH:6][CH:7]=[CH:8][CH:9]=[C:10]12>>[ClH:1].[N:2]1([C@H:11]([C@@H:12]([CH2:13][NH:14][CH3:15])[OH:16])[C:17]2=[CH:22][CH:21]=[C:20]([O:23][C:24]([F:25])([F:26])[F:27])[CH:19]=[CH:18]2)[CH:3]=[CH:4][C:5]2=[CH:6][CH:7]=[CH:8][CH:9]=[C:10]12"
    smi = "CC(C)(C)OC(=O)[N:21]1[CH2:20][CH:19]([NH:18][C:16]([C:12]2=[CH:13][CH:14]=[C:15]3[N:7]([CH2:6][C:5]4=[CH:4][CH:3]=[C:2]([F:1])[CH:34]=[CH:33]4)[C:8]([CH3:32])=[C:9]([CH3:31])[C:10]3=[CH:11]2)=[O:17])[CH2:23][CH2:22]1.[ClH:35]>>[ClH:35].[F:1][C:2]1=[CH:34][CH:33]=[C:5]([CH2:6][N:7]2[C:8]([CH3:32])=[C:9]([CH3:31])[C:10]3=[CH:11][C:12]([C:16](=[O:17])[NH:18][CH:19]4[CH2:20][NH:21][CH2:22][CH2:23]4)=[CH:13][CH:14]=[C:15]23)[CH:4]=[CH:3]1"
    print(TemplateExtractor.rxn_smiles_to_rxn_temp_smarts(smi, RxnTemplateType.EXTENDED, True))
