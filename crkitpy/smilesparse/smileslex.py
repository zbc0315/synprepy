# smileslex.py. This file automatically created by PLY (version 3.11). Don't edit!
_tabversion   = '3.10'
_lextokens    = set(('ATOM', 'ATOM_H', 'BOND', 'CHARGE', 'COLON', 'CYCLE_NUM', 'DOT', 'LPAREN', 'LSQUAR', 'MAP_NUM', 'NUM', 'RPAREN', 'RSQUAR'))
_lexreflags   = 64
_lexliterals  = ''
_lexstateinfo = {'INITIAL': 'inclusive', 'inquare': 'inclusive'}
_lexstatere   = {'INITIAL': [('(?P<t_begin_inquare>\\[)|(?P<t_ATOM>Cl|Br|B|C|N|O|F|P|S|I|b|c|n|o|p|s|\\*)|(?P<t_CHARGE>\\+\\d+|\\-\\d+|\\++|\\-+)|(?P<t_CYCLE_NUM>[0-9]+\\%)|(?P<t_NUM>[0-9]+)|(?P<t_BOND>[-=#:])|(?P<t_DOT>\\.)|(?P<t_RPAREN>\\))|(?P<t_LPAREN>\\()|(?P<t_COLON>:)', [None, ('t_begin_inquare', 'begin_inquare'), ('t_ATOM', 'ATOM'), ('t_CHARGE', 'CHARGE'), ('t_CYCLE_NUM', 'CYCLE_NUM'), ('t_NUM', 'NUM'), (None, 'BOND'), (None, 'DOT'), (None, 'RPAREN'), (None, 'LPAREN'), (None, 'COLON')])], 'inquare': [('(?P<t_inquare_end>\\])|(?P<t_inquare_ATOM>He|Li|Be|Ne|Na|Mg|Al|Si|Ar|K|Ca|Sc|Ti|V|Cr|Mn|Fe|Co|Ni|Cu|Zn|Ga|Ge|As|Se|Kr|Rb|Sr|Y|Zr|Nb|Mo|Tc|Ru|Rh|Pd|Ag|Cd|In|Sn|Sb|Te|Xe|Cs|Ba|La|Ce|Pr|Nd|Pm|Sm|Eu|Gd|Tb|Dy|Ho|Er|Tm|Yb|Lu|Hf|Ta|W|Re|Os|Ir|Pt|Au|Hg|Tl|Pb|Bi|Po|At|Rn|Fr|Ra|Ac|Th|Pa|U|Np|Pu|Am|Cm|Bk|Cf|Es|Fm|Md|No|Lr|Rf|Db|Sg|Bh|Hs|Mt|Ds|Rg|Cn|Nh|Fl|Mc|Lv|Ts|Og|Uun|Uuu|Uub|Uut|Uuq|Uup|Uuh|Uus|Uuo|si|as|se|te|ARO)|(?P<t_inquare_ATOM_H>H)|(?P<t_inquare_MAP_NUM>\\:[0-9]+)', [None, ('t_inquare_end', 'end'), ('t_inquare_ATOM', 'ATOM'), ('t_inquare_ATOM_H', 'ATOM_H'), ('t_inquare_MAP_NUM', 'MAP_NUM')]), ('(?P<t_begin_inquare>\\[)|(?P<t_ATOM>Cl|Br|B|C|N|O|F|P|S|I|b|c|n|o|p|s|\\*)|(?P<t_CHARGE>\\+\\d+|\\-\\d+|\\++|\\-+)|(?P<t_CYCLE_NUM>[0-9]+\\%)|(?P<t_NUM>[0-9]+)|(?P<t_BOND>[-=#:])|(?P<t_DOT>\\.)|(?P<t_RPAREN>\\))|(?P<t_LPAREN>\\()|(?P<t_COLON>:)', [None, ('t_begin_inquare', 'begin_inquare'), ('t_ATOM', 'ATOM'), ('t_CHARGE', 'CHARGE'), ('t_CYCLE_NUM', 'CYCLE_NUM'), ('t_NUM', 'NUM'), (None, 'BOND'), (None, 'DOT'), (None, 'RPAREN'), (None, 'LPAREN'), (None, 'COLON')])]}
_lexstateignore = {'INITIAL': '', 'inquare': ''}
_lexstateerrorf = {'INITIAL': 't_error', 'inquare': 't_error'}
_lexstateeoff = {}