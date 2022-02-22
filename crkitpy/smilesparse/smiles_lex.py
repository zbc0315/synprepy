# -*- coding: utf-8 -*-
# @Time     : 2020/5/11 16:19
# @Author   : zbc
# @Email    : zhangbc0315@outlook.com
# @FileName : smiles_lex.py
# @Software : PyCharm

import ply.lex as lex

from crkitpy.molgraph.atom import Atom

states = [
    ('inquare', 'inclusive'),
]

tokens = [
    'DOT',
    'LPAREN',
    'RPAREN',
    'LSQUAR',
    'RSQUAR',
    'CHARGE',
    'MAP_NUM',
    'CYCLE_NUM',
    'NUM',
    'BOND',
    'ATOM',
    'ATOM_H',
    'COLON',
]


def SmilesLex():

    # t_AT = r'@'
    t_LPAREN = r'\('
    t_RPAREN = r'\)'
    t_DOT = r'\.'
    t_COLON = r':'

    t_BOND = r'[-=#:]'

    r_atom = r'Cl|Br|B|C|N|O|F|P|S|I|b|c|n|o|p|s|\*'
    r_inquare_atom = r'He|Li|Be|Ne|Na|Mg|Al|Si|Ar|' \
                     r'K|Ca|Sc|Ti|V|Cr|Mn|Fe|Co|Ni|Cu|Zn|Ga|Ge|As|Se|Kr|' \
                     r'Rb|Sr|Y|Zr|Nb|Mo|Tc|Ru|Rh|Pd|Ag|Cd|In|Sn|Sb|Te|Xe|' \
                     r'Cs|Ba|' \
                     r'La|Ce|Pr|Nd|Pm|Sm|Eu|Gd|Tb|Dy|Ho|Er|Tm|Yb|Lu|' \
                     r'Hf|Ta|W|Re|Os|Ir|Pt|Au|Hg|Tl|Pb|Bi|Po|At|Rn|' \
                     r'Fr|Ra|' \
                     r'Ac|Th|Pa|U|Np|Pu|Am|Cm|Bk|Cf|Es|Fm|Md|No|Lr|' \
                     r'Rf|Db|Sg|Bh|Hs|Mt|Ds|Rg|Cn|Nh|Fl|Mc|Lv|Ts|Og|Uun|Uuu|Uub|Uut|Uuq|Uup|Uuh|Uus|Uuo|' \
                     r'si|as|se|te|ARO'
    r_inquare_atom_h = r'H'

    def t_begin_inquare(t):
        r"""\["""
        t.lexer.begin('inquare')
        t.type = 'LSQUAR'
        return t

    def t_inquare_end(t):
        r"""\]"""
        t.lexer.begin('INITIAL')
        t.type = 'RSQUAR'
        return t

    @lex.TOKEN(r_atom)
    def t_ATOM(t):
        atom = Atom.instance_by_symbol(t.value)
        if t.lexer.lexstate == 'inquare':
            atom.implicit_hs_num = 0
        init_aromatic(atom)
        t.value = atom
        return t

    @lex.TOKEN(r_inquare_atom)
    def t_inquare_ATOM(t):
        atom = Atom.instance_by_symbol(t.value)
        atom.implicit_hs_num = 0
        init_aromatic(atom)
        t.value = atom
        return t

    @lex.TOKEN(r_inquare_atom_h)
    def t_inquare_ATOM_H(t):
        atom = Atom.instance_by_symbol(t.value)
        t.value = atom
        return t

    def t_CHARGE(t):
        r"""\+\d+|\-\d+|\++|\-+"""

        if len(t.value) == 1:
            t.value = 1 if t.value == '+' else -1
        elif t.value[-1].isdigit():
            t.value = int(t.value)
        else:
            t.value = len(t.value) if t.value[0] == '+' else -len(t.value)
        return t

    def t_CYCLE_NUM(t):
        r"""[0-9]+\%"""
        t.value = int(t.value[:-1])
        return t

    def t_NUM(t):
        r"""[0-9]+"""
        t.value = int(t.value)
        return t

    def t_inquare_MAP_NUM(t):
        r"""\:[0-9]+"""
        t.value = int(t.value[1:])
        return t

    def t_error(t):
        print('Syntax error is input in lex!')

    def init_aromatic(atom: Atom) -> None:
        if atom.get_symbol().islower():
            atom.is_aromatic = True
            atom.set_symbol(atom.get_symbol().capitalize())

    return lex.lex(optimize=1, lextab='smileslex')


if __name__ == '__main__':
    m = SmilesLex()
    m.input("CccccCH([Sc][Cl][PH3])[Mn]ClScCC:3C2")
    # m.input("BCCl[He:2](CCC)[B][Be]")
    while True:
        tok = m.token()
        if tok is None:
            break
        print(tok)
