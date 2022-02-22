#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Time    : 2021/9/13 13:59
# @Author  : zhangbc0315@outlook.com
# @File    : rxn_utils.py
# @Software: PyCharm

import pandas as pd
from rdkit.Chem import AllChem

from mcts_lite.mol_utils import MolUtils
from mcts_lite.data.db_spectators import DbSpectators


class RxnUtils:

    rxn_df = None

    @classmethod
    def save(cls, fp: str):
        if cls.rxn_df is not None:
            cls.rxn_df.to_csv(fp, sep='\t', encoding='utf-8')

    @classmethod
    def init(cls):
        cls.rxn_df = pd.DataFrame(columns=['p_mid', 'rs_mids', 'cats_oh_idxes', 'sols_oh_idxes', 'prob', 'tid'])

    @classmethod
    def query_rxns_by_p_mid(cls, p_mid):
        res = cls.rxn_df.query(f'p_mid=={p_mid}')
        return list(res.index)

    @classmethod
    def query_by_rid(cls, rid):
        return cls.rxn_df.loc[rid, ]

    @classmethod
    def add_rxn(cls, p_mid, rs_mids, cats_oh_idxes, sols_oh_idxes, prob, tid):
        cls.rxn_df = cls.rxn_df.append({'p_mid': p_mid, 'rs_mids': set(rs_mids), 'cats_oh_idxes': cats_oh_idxes,
                                        'sols_oh_idxes': sols_oh_idxes, 'prob': prob, 'tid': tid}, ignore_index=True)
        return len(cls.rxn_df) - 1

    @classmethod
    def get_rxn_dict(cls, rid):
        rxn = cls.query_by_rid(rid)
        res = {'p_inchi': MolUtils.query_by_mid(rxn.p_mid).inchi,
               'rs_inchis': [MolUtils.query_by_mid(mid).inchi for mid in rxn.rs_mids],
               'cats_codes': DbSpectators.cats_idxes_to_codes(map(lambda x: x-1, rxn.cats_oh_idxes)),
               'sols_codes': DbSpectators.sols_idxes_to_codes(map(lambda x: x-1, rxn.sols_oh_idxes))}
        return res

    @classmethod
    def get_rxn_smiles(cls, rid):
        rxn_data = cls.query_by_rid(rid)
        rs_inchis = [MolUtils.query_by_mid(mid).inchi for mid in rxn_data.rs_mids]
        p_inchi = MolUtils.query_by_mid(rxn_data.p_mid).inchi
        rxn = AllChem.ChemicalReaction()
        rxn.AddProductTemplate(AllChem.MolFromInchi(p_inchi))
        for r_inchi in rs_inchis:
            rxn.AddReactantTemplate(AllChem.MolFromInchi(r_inchi))
        rxn_smiles = AllChem.ReactionToSmiles(rxn)
        cats = cls._list_to_str(map(lambda x: x-1, rxn_data.cats_oh_idxes))
        sols = cls._list_to_str(map(lambda x: x-1, rxn_data.sols_oh_idxes))
        return f'{rxn_smiles} |c| {cats} |s| {sols} |t| {rxn_data.tid}'

    @classmethod
    def _list_to_str(cls, li):
        l_str = [str(n) for n in li]
        return '.'.join(l_str)


if __name__ == "__main__":
    RxnUtils.init()
    RxnUtils.add_rxn(1, {2, 3}, [], [], 0.5)
    RxnUtils.add_rxn(2, [2, 3], [], [], 0.5)
    RxnUtils.add_rxn(1, [3], [], [], 0.5)
    print(RxnUtils.query_rxns_by_p_mid(1))
