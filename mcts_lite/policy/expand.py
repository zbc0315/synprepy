#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Time    : 2021/9/1 16:13
# @Author  : zhangbc0315@outlook.com
# @File    : expand.py
# @Software: PyCharm

from rdkit.Chem import AllChem

from mcts_lite.policy.basic_temp_selector import BasicTempSelector
from mcts_lite.data.db_rxn_template import DbRxnTemplate
from mcts_lite.policy.check_reactants import CheckReactants
from mcts_lite.policy.sane_determiner import SaneDeterminer
from mcts_lite.node_utils import NodeUtils
from mcts_lite.rxn_utils import RxnUtils
from mcts_lite.mol_utils import MolUtils
from crkitpy.smatesparse.smates_parser import SmatesParser
from crkitpy.smilesparse.rdkit_smiles_parser import RdkitSmilesParser
from crkitpy.rxnhandler.rxn_template_apply import RxnTemplateApply
from log_utils.logger import Logging


class Expand:

    db_rxn_template = DbRxnTemplate()
    cat_num = 325
    sol_num = 321

    @classmethod
    def expand_by_nid(cls, nid, target_inchi):
        inchi, p_mid = NodeUtils.get_unsolved_inchi_and_mid(nid)
        existed_rids = RxnUtils.query_rxns_by_p_mid(p_mid)
        intermediate_products = NodeUtils.get_intermediate_products(nid, p_mid)
        if len(existed_rids) == 0:
            try:
                rids = cls.expand_by_inchi(inchi, target_inchi, p_mid, intermediate_products)
            except Exception as e:
                raise e
                Logging.log.warning(f'expand error: {e}')
                NodeUtils.set_no_expand(nid, True)
                NodeUtils.set_no_select(nid, True)
                return
            if len(rids) == 0:
                NodeUtils.set_no_expand(nid, True)
                NodeUtils.set_no_select(nid, True)
            NodeUtils.create_child_nodes_by_rids(nid, rids)
        else:
            NodeUtils.create_child_nodes_by_rids(nid, existed_rids)

    @classmethod
    def expand_by_inchi(cls, inchi, target_inchi, p_mid, intermediate_products):
        Logging.log.info(f'Expand - target inchi: {inchi}')
        rd_product = AllChem.MolFromInchi(inchi)

        if rd_product.GetNumAtoms() <= 1:
            return None

        total_prob = 0
        probs = []
        yielded_reactants_inchis = []
        rd_reactants_list = []
        tids = []
        smateses = []
        num_right_tid = 0
        for n, (tid, prob) in enumerate(BasicTempSelector.predict(rd_product)):
            # if inchi == 'InChI=1S/C18H21ClN5O7P/c1-5-29-32(27,30-6-2)31-10-23-8-12(16(25)17(19)26)14-13(28-4)7-20-18(15(14)23)24-9-21-11(3)22-24/h7-9H,5-6,10H2,1-4H3' and tid == 26700:
            #     print(1)
            have_reactants = False
            right = 0
            wrong = 0
            smates, cover_score = cls.db_rxn_template.search_basic_smates_by_tid(tid)
            if smates is None:
                Logging.log.warning(f'Can not find smates for tid: {tid}')
                continue
            num_right_tid += 1
            if num_right_tid >= 10:
                break
            Logging.log.info(f'Expand - {tid} - {smates}')
            product = RdkitSmilesParser.rdmol_to_mol(rd_product)
            for reactants in RxnTemplateApply.product_to_reactants(smates, product):
                reactants_inchis, rd_reactants, is_right = CheckReactants.check(reactants, inchi, target_inchi,
                                                                                yielded_reactants_inchis,
                                                                                intermediate_products)
                if not is_right:
                    wrong += 1
                    continue
                have_reactants = True
                right += 1
                yielded_reactants_inchis.append(set(reactants_inchis))
                rd_reactants_list.append(rd_reactants)
                tids.append(tid)
                smateses.append(smates)
                probs.append(min(prob + cover_score, 1))

            total_prob += prob
            if total_prob >= 0.995:
                Logging.log.debug(f'总概率达到: {total_prob}, 提前结束expand')
                break

            if wrong > 0:
                Logging.log.debug(f'Expand - tid: {tid} - {right}/{wrong} - smates: {smates}')
        if len(rd_reactants_list) == 0:
            return []
        return list(cls._get_mcts_rxns(rd_reactants_list, rd_product, tids, smateses, yielded_reactants_inchis, p_mid, probs))

    @classmethod
    def _get_mcts_rxns(cls, rd_reactants_list, rd_product, tids, smateses, yielded_reactants_inchis, p_mid, probs):
        for n, cat_oh_idxes, sol_oh_idxes, tid in cls.sane_reactions(rd_reactants_list, rd_product, tids, smateses):
            reactants_inchis = yielded_reactants_inchis[n]
            reactants_mids = list(MolUtils.add_or_query_inchis(reactants_inchis))
            rid = RxnUtils.add_rxn(p_mid, reactants_mids, cat_oh_idxes, sol_oh_idxes, probs[n], tid)
            yield rid

    # ============================
    # --    Is Reaction Sane    --

    @classmethod
    def _sort_idxes(cls, oh_idxes, sorted_all_idxes):
        res = []
        for idx in sorted_all_idxes:
            if idx+1 in oh_idxes:
                res.append(idx+1)
                oh_idxes.remove(idx+1)
        for idx in oh_idxes:
            res.append(idx)
        return res

    @classmethod
    def _simplify_cat_sol(cls, cat_sol_oh_idx_list, sorted_cat_idxes, sorted_sol_idxes):
        cat_oh_idxes = set()
        sol_oh_idxes = set()
        for cat, sol in cat_sol_oh_idx_list:
            cat_oh_idxes.add(cat)
            sol_oh_idxes.add(sol)
        cat_oh_idxes = cls._sort_idxes(cat_oh_idxes, sorted_cat_idxes)
        sol_oh_idxes = cls._sort_idxes(sol_oh_idxes, sorted_sol_idxes)

        return cat_oh_idxes, sol_oh_idxes

    @classmethod
    def sane_reactions(cls, rd_reactants_list, rd_product, tids, smateses, consider_cat_sol: bool = True):
        if consider_cat_sol:
            for n, (is_sane, cat_sol_oh_idx_list, tid, sorted_cat_idxes, sorted_sol_idxes) in \
                    enumerate(SaneDeterminer.is_sane_consider_spec(rd_reactants_list,
                                                                   rd_product,
                                                                   tids,
                                                                   smateses)):
                if not is_sane:
                    continue
                cat_oh_idxes, sol_oh_idxes = cls._simplify_cat_sol(cat_sol_oh_idx_list, sorted_cat_idxes, sorted_sol_idxes)
                yield n, cat_oh_idxes, sol_oh_idxes, tid


if __name__ == "__main__":
    pass
