# -*- coding: utf-8 -*-
# @Time     : 2020/6/4 15:10
# @Author   : zbc
# @Email    : zhangbc0315@outlook.com
# @FileName : smates_parser.py
# @Software : PyCharm

from crkitpy.molgraph.reaction import Reaction
from crkitpy.molgraph.molecule import Molecule, Atom
from crkitpy.smilesparse.smiles_parser import SmilesParser
from crkitpy.const import ATOM_PROP
from crkitpy.chemdata.atomic_data import AtomicData


class SmatesParser:

    @classmethod
    def template_from_smates(cls, smates: str) -> Reaction:
        if ' | ' in smates:
            rxn_info, other_info = smates.split(' | ')
        else:
            rxn_info = smates
            other_info = None
        reactant_smates, product_smates = rxn_info.split('>>')
        rxn_template = Reaction()
        rxn_template.add_reactant(SmilesParser.parse_to_mol_by_lex(reactant_smates))
        rxn_template.add_product(SmilesParser.parse_to_mol_by_lex(product_smates))
        cls._detect_required_and_max_valence(rxn_template)
        if other_info is not None:
            cls._extract_rxn_other_info(rxn_template, other_info)
        return rxn_template

    @classmethod
    def _extract_rxn_other_info(cls, rxn: Reaction, other_info: str):
        """
        :param rxn:
        :param other_info: f'PM{map_num}:{min_hs},{max_hs};'

        """
        mols = {'R': rxn.get_reactants()[0],
                'P': rxn.get_products()[0]}
        for atom_other_info in other_info.split(';'):
            mol = mols[atom_other_info[0]]
            cls._extract_atom_other_info(mol, atom_other_info)

    @classmethod
    def _extract_atom_other_info(cls, mol: Molecule, atom_other_info: str):
        atom_key, atom_values = atom_other_info.split(':')
        atoms = cls._get_atoms(mol, atom_key)
        min_hs, max_hs = atom_values.split(',')
        if len(min_hs) > 0:
            [atom.set_property(ATOM_PROP.K_MIN_HS, int(min_hs)) for atom in atoms]
        if len(max_hs) > 0:
            [atom.set_property(ATOM_PROP.K_MAX_HS, int(max_hs)) for atom in atoms]

    @classmethod
    def _get_atoms(cls, mol: Molecule, atom_key: str) -> [Atom]:
        if atom_key[1] == 'M':
            map_num = int(atom_key[2:])
            return list(mol.get_atoms_by_map_num(map_num))
        if atom_key[1] == 'I':
            raise ValueError("Haven't Set Atom Index")

    @classmethod
    def _detect_required_and_max_valence(cls, rxn_template: Reaction):
        """
        for p_atom in rxn_template_selector.get_products()[0].get_atoms():
            r_atom = rxn_template_selector.get_reactants()[0].get_one_atom_by_map_num(p_atom.map_num)
            p_require_valence = r_atom.get_bonds_weight_without_h() - p_atom.get_bonds_weight_without_h()
            p_atom.set_property(ATOM_PROP.K_REQUIRE_VALENCE, p_require_valence)
            r_atom.set_property(ATOM_PROP.K_REQUIRE_VALENCE, -p_require_valence)

            max_valence = max(AtomicData.get_valences(p_atom.get_symbol()))
            if max_valence < 0:
                max_valence = 0
            elif max_valence < r_atom.get_valence():
                max_valence = r_atom.get_valence()
            p_atom.set_property(ATOM_PROP.K_MAX_VALENCE, max_valence)
            r_atom.set_property(ATOM_PROP.K_MAX_VALENCE, max_valence)
        """
        for r_atom in rxn_template.get_reactants()[0].get_atoms():
            if r_atom.is_r_group():
                continue
            max_valence = max(AtomicData.get_valences(r_atom.get_symbol()))
            if max_valence < 0:
                max_valence = 0
            elif max_valence < r_atom.get_valence():
                max_valence = r_atom.get_valence()

            p_atom = rxn_template.get_products()[0].get_one_atom_by_map_num(r_atom.map_num)
            if p_atom is None:
                r_atom.set_property(ATOM_PROP.K_REQUIRE_VALENCE, 0)
                r_atom.set_property(ATOM_PROP.K_MAX_VALENCE, max_valence)
            else:
                p_require_valence = r_atom.get_bonds_weight_without_h() - p_atom.get_bonds_weight_without_h()
                p_atom.set_property(ATOM_PROP.K_REQUIRE_VALENCE, p_require_valence)
                r_atom.set_property(ATOM_PROP.K_REQUIRE_VALENCE, -p_require_valence)

                p_atom.set_property(ATOM_PROP.K_MAX_VALENCE, max_valence)
                r_atom.set_property(ATOM_PROP.K_MAX_VALENCE, max_valence)

