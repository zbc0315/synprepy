# -*- coding: utf-8 -*-
# @Time     : 2020/5/13 14:21
# @Author   : zbc
# @Email    : zhangbc0315@outlook.com
# @FileName : smiles_yacc.py
# @Software : PyCharm

import ply.yacc as yacc

from crkitpy.smilesparse.smiles_lex import SmilesLex, tokens
from crkitpy.molgraph.molecule import Molecule, Atom
from crkitpy.const import BOND_TYPE, ATOM_PROP


def SmilesYacc():

    # 起始是原子
    def p_atom(p):
        """expression : ATOM"""
        molecule = Molecule()
        atom = p[1]
        molecule.add_atom(atom)
        p[0] = molecule
        set_activate_atom(atom)

    # 起始是'['，其后紧跟原子
    def p_atom_in_square(p):
        """expression : LSQUAR ATOM
                      | LSQUAR ATOM_H"""
        molecule = Molecule()
        atom = p[2]
        molecule.add_atom(atom)
        p[0] = molecule
        set_activate_atom(atom)

    def p_bond_atom_h_start_by_lsquar(p):
        """expression : expression LSQUAR ATOM_H"""
        nonlocal activate_link_atom
        nonlocal activate_prop_atom
        p[0] = p[1]
        atom_h = p[3]
        p[0].add_atom(atom_h)
        if activate_link_atom is not None:
            activate_link_atom.link_bond_by_type(current_bond_type, atom_h)
        else:
            activate_link_atom = atom_h
        activate_prop_atom = atom_h

    def p_dot(p):
        """expression : expression DOT"""
        nonlocal activate_link_atom
        nonlocal activate_prop_atom
        activate_link_atom = None
        activate_prop_atom = None
        p[0] = p[1]

    def p_charge(p):
        """expression : expression CHARGE"""
        nonlocal activate_prop_atom
        p[0] = p[1]
        activate_prop_atom.charge = p[2]

    # 化学键
    def p_bond(p):
        """expression : expression BOND"""
        nonlocal current_bond_type
        p[0] = p[1]
        if p[2] == '-':
            current_bond_type = BOND_TYPE.SINGLE
        elif p[2] == '=':
            current_bond_type = BOND_TYPE.DOUBLE
        elif p[2] == '#':
            current_bond_type = BOND_TYPE.TRIPLE
        elif p[2] == ':':
            current_bond_type = BOND_TYPE.AROMATIC

    def p_bond_atom(p):
        """expression : expression ATOM"""
        nonlocal activate_link_atom
        nonlocal current_bond_type
        p[0] = p[1]
        atom = p[2]
        p[0].add_atom(atom)
        if activate_link_atom is not None:
            link_bond_(activate_link_atom, atom)
            # if atom.is_aromatic and activate_link_atom.is_aromatic:
            #     activate_link_atom.link_bond_by_type(BOND_TYPE.AROMATIC, atom)
            # else:
            #     activate_link_atom.link_bond_by_type(current_bond_type, atom)
        set_activate_atom(atom)
        current_bond_type = BOND_TYPE.SINGLE

    def p_bond_many_atom_h(p):
        """expression : expression ATOM_H NUM"""
        nonlocal activate_link_atom
        p[0] = p[1]
        activate_link_atom.explicit_hs_num = p[3]
        activate_link_atom.set_implicit_valence(p[3])

    def p_bond_one_atom_h(p):
        """expression : expression ATOM_H"""
        nonlocal activate_link_atom
        p[0] = p[1]
        activate_link_atom.explicit_hs_num = 1
        activate_link_atom.set_implicit_valence(1)

    def p_begin_paren(p):
        """expression : expression LPAREN"""
        nonlocal branch_atoms
        nonlocal activate_link_atom
        branch_atoms.append(activate_link_atom)
        p[0] = p[1]

    def p_end_paren(p):
        """expression : expression RPAREN"""
        nonlocal branch_atoms
        set_activate_atom(branch_atoms[-1])
        branch_atoms.pop(-1)
        p[0] = p[1]

    def p_map_num(p):
        """expression : expression MAP_NUM"""
        nonlocal activate_prop_atom
        p[0] = p[1]
        activate_prop_atom.map_num = p[2]

    def p_complex_cycle_num(p):
        """expression : expression CYCLE_NUM"""
        nonlocal activate_link_atom
        nonlocal cycle_nums_atoms
        p[0] = p[1]
        cycle_num = p[2]
        if cycle_num not in cycle_nums_atoms.keys():
            cycle_nums_atoms[cycle_num] = activate_link_atom
        else:
            link_bond_(activate_link_atom, cycle_nums_atoms[cycle_num])
            # activate_link_atom.link_bond_by_type(current_bond_type, cycle_nums_atoms[cycle_num])
            cycle_nums_atoms.pop(cycle_num)

    def p_simple_cycle_num(p):
        """expression : expression NUM"""
        nonlocal activate_link_atom
        nonlocal cycle_nums_atoms
        p[0] = p[1]
        cycle_nums = [int(n) for n in str(p[2])]
        for cycle_num in cycle_nums:
            if cycle_num not in cycle_nums_atoms.keys():
                cycle_nums_atoms[cycle_num] = activate_link_atom
            else:
                link_bond_(activate_link_atom, cycle_nums_atoms[cycle_num])
                # activate_link_atom.link_bond_by_type(current_bond_type, cycle_nums_atoms[cycle_num])
                cycle_nums_atoms.pop(cycle_num)

    def p_only_square(p):
        """expression : expression RSQUAR
                      | expression LSQUAR"""
        p[0] = p[1]

    def p_colon(p):
        """expression : expression COLON"""
        p[0] = p[1]

    def p_error(p):
        print(f'Syntax error in input in yacc: {p}!')

    def set_activate_atom(atom: Atom):
        nonlocal activate_link_atom
        nonlocal activate_prop_atom
        activate_link_atom = atom
        activate_prop_atom = atom

    def link_bond_(begin_atom: Atom, end_atom: Atom):
        nonlocal current_bond_type
        # if begin_atom.is_aromatic and end_atom.is_aromatic:
        #     begin_atom.link_bond_by_type(BOND_TYPE.AROMATIC, end_atom)
        # else:
        begin_atom.link_bond_by_type(current_bond_type, end_atom)

    activate_link_atom: Atom = None
    activate_prop_atom: Atom = None
    current_bond_type = BOND_TYPE.SINGLE
    branch_atoms = []
    cycle_nums_atoms = {}
    return yacc.yacc()


if __name__ == '__main__':
    smi = "[CH3:9][CH:8]=[C:5]12([CH3:6])[C:4]2%3%(=[CH2:7])[CH:3]([C:10]#[CH:11])[CH:2]([CH3:1])[CH3:12]C:C"
    smiles_lex = SmilesLex()
    smiles_yacc = SmilesYacc()
    result = smiles_yacc.parse(smiles_lex.input(smi))
    result.draw()
