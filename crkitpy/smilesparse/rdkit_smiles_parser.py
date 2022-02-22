# -*- coding: utf-8 -*-
# @Time     : 2020/4/21 10:13
# @Author   : zbc
# @Email    : zhangbc0315@outlook.com
# @FileName : rdkit_smiles_parser.py
# @Software : PyCharm

from rdkit.Chem import AllChem
from rdkit.Chem import rdChemReactions

from crkitpy.molgraph.reaction import Reaction
from crkitpy.molgraph.molecule import Molecule
from crkitpy.molgraph.atom import Atom
from crkitpy.molgraph.bond import Bond
from crkitpy.errors.smiles_error import RdkitSmilesError


class RdkitSmilesParser:

    @classmethod
    def parse_to_rxn(cls, smiles: str) -> Reaction:
        split_rxn_smiles = smiles.split('>')
        reactant_smiles_list = split_rxn_smiles[0].split('.')
        product_smiles_list = split_rxn_smiles[-1].split('.')
        reaction = Reaction()
        for reactant_smiles in reactant_smiles_list:
            reaction.add_reactant(cls.parse_to_mol(reactant_smiles))
        for product_smiles in product_smiles_list:
            reaction.add_product(cls.parse_to_mol(product_smiles))
        return reaction

    @classmethod
    def parse_to_mol(cls, smiles: str) -> Molecule:
        rdmol = AllChem.MolFromSmiles(smiles)
        if rdmol is None:
            raise RdkitSmilesError(f'rdkit can not parse this smiles: {smiles}')
        return cls.rdmol_to_mol(rdmol)

    @classmethod
    def rdmol_to_mol(cls, rdmol) -> Molecule:
        molecule = Molecule()
        atoms_map = {}
        for rdatom in rdmol.GetAtoms():
            atom = cls._rdatom_to_atom(rdatom)
            atoms_map[rdatom.GetIdx()] = atom
            molecule.add_atom(atom)
        for rdbond in rdmol.GetBonds():
            begin_atom = atoms_map[rdbond.GetBeginAtomIdx()]
            end_atom = atoms_map[rdbond.GetEndAtomIdx()]
            molecule.add_bond(cls._rdbond_to_bond(rdbond), begin_atom, end_atom)
        return molecule

    @staticmethod
    def _rdatom_to_atom(rdatom):
        atom = Atom.instance_by_symbol(rdatom.GetSymbol())
        atom.set_id(rdatom.GetIdx())  # TODO 立刻取消！！！
        atom.map_num = rdatom.GetAtomMapNum()
        atom.degree = rdatom.GetDegree()
        atom.set_atomic_num(rdatom.GetAtomicNum())
        atom.is_aromatic = rdatom.GetIsAromatic()
        atom.charge = rdatom.GetFormalCharge()
        atom.set_explicit_valence(rdatom.GetExplicitValence())
        atom.set_implicit_valence(rdatom.GetImplicitValence())
        atom.explicit_hs_num = rdatom.GetNumExplicitHs()
        atom.implicit_hs_num = rdatom.GetNumImplicitHs()
        return atom

    @staticmethod
    def _rdbond_to_bond(rdbond):
        bond = Bond.instance_by_bond_type(str(rdbond.GetBondType()))
        return bond
