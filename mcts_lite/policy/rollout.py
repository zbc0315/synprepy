#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Time    : 2021/9/7 16:29
# @Author  : zhangbc0315@outlook.com
# @File    : rollout.py
# @Software: PyCharm

from mcts_lite.policy.comp_temp_selector import CompTempSelector
from mcts_lite.data.db_rxn_template import DbRxnTemplate
from mcts_lite.policy.check_reactants import CheckReactants
from crkitpy.smatesparse.smates_parser import SmatesParser
from crkitpy.inchiparse.inchi_parser import InchiParser
from crkitpy.rxnhandler.rxn_template_apply import RxnTemplateApply
from mcts_lite.data.db_cac import DbCac
from mcts_lite.mol_utils import MolUtils
from mcts_lite.node_utils import NodeUtils
from log_utils.logger import Logging


class Rollout:

    @classmethod
    def _get_cac_score(cls, inchis):
        cac = []
        no_cac = []
        for inchi in inchis:
            if DbCac.is_cac(inchi):
                cac.append(inchi)
            else:
                no_cac.append(inchi)
        return len(cac) / (len(cac) + len(no_cac)), cac, no_cac

    @classmethod
    def _rollout_for_mol(cls, inchi, comp_smates):
        """ 模拟一个分子的未来发展，尝试它的每一个模板，尝试每一种原料集合，
        得分 = 商业可获得原料数 / 原料数，
        返回得分最高的原料，
        如果每一种原料集合都不合法，则返回分数为0，返回的原料数目为0
        """
        Logging.log.debug(f'Rollout - Target inchi - {inchi}')
        best_score = None
        cac_inchis = []
        no_cac_inchis = []
        for smates in comp_smates:
            right = 0
            wrong = 0
            product = InchiParser.mol_from_inchi(inchi)
            for reactants in RxnTemplateApply.product_to_reactants(smates, product):
                reactant_inchis, _, is_right = CheckReactants.check(reactants, inchi, inchi, [], [])
                if not is_right:
                    wrong += 1
                    continue
                right += 1
                score, cac, no_cac = cls._get_cac_score(reactant_inchis)
                if best_score is None or score >= best_score:
                    cac_inchis = cac
                    no_cac_inchis = no_cac
                    best_score = score
                if best_score == 1:
                    break
            if wrong > 0:
                Logging.log.debug(f'Rollout - Smates - {smates}')
            if best_score == 1:
                break
        if best_score is None:
            return 0, [], [inchi]
        return best_score, cac_inchis, no_cac_inchis

    @classmethod
    def rollout_by_inchis(cls, all_cac_inchis, all_no_cac_inchis, rest_step_num):
        # all_cac = set()
        new_all_no_cac_inchis = set()
        if len(all_no_cac_inchis) == 0:
            return all_cac_inchis, set(all_no_cac_inchis)
        for n, comp_tids in enumerate(CompTempSelector.batch_predict_by_inchis(all_no_cac_inchis)):
            inchi = all_no_cac_inchis[n]
            comp_smates = [DbRxnTemplate.search_comp_smates_by_tid(tid) for tid in comp_tids]
            score, cac_inchis, no_cac_inchis = cls._rollout_for_mol(inchi, comp_smates)
            all_cac_inchis = all_cac_inchis.union(set(cac_inchis))
            new_all_no_cac_inchis = new_all_no_cac_inchis.union(set(no_cac_inchis))
            # all_cac = all_cac.union(set(cac_inchis))
            # all_no_cac = all_no_cac.union(set(no_cac_inchis))
            rest_step_num -= 1
            if rest_step_num == 0:
                new_all_no_cac_inchis = new_all_no_cac_inchis.union(set(all_no_cac_inchis[n+1:]))
                return all_cac_inchis, new_all_no_cac_inchis
        if rest_step_num == 0:
            return all_cac_inchis, new_all_no_cac_inchis
        else:
            return cls.rollout_by_inchis(all_cac_inchis, list(new_all_no_cac_inchis), rest_step_num)
        # return all_cac, all_no_cac

    @classmethod
    def rollout_by_mols(cls, mids, rest_step_num):
        inchis = []
        all_cac = set()
        for mid in mids:
            mol_data = MolUtils.query_by_mid(mid)
            if not mol_data.is_cac:
                inchis.append(mol_data.inchi)
            else:
                all_cac.add(mol_data.inchi)
        rest_step_num = min(rest_step_num, len(inchis))
        all_cac, no_cac = cls.rollout_by_inchis(all_cac, inchis, rest_step_num)
        # all_cac = all_cac.union(cac)
        all_no_cac = no_cac
        if len(all_cac) + len(all_no_cac) == 0:
            return 0
        return len(all_cac) / (len(all_cac) + len(all_no_cac))

    @classmethod
    def rollout_by_node(cls, nid):
        rest_step = NodeUtils.get_rest_step(nid)
        return cls.rollout_by_mols(NodeUtils.query_by_nid(nid).mids, rest_step)


def test():
    inchis = ['InChI=1S/C17H28/c1-2-9-15-13(5-1)11-14-8-3-6-12-7-4-10-16(15)17(12)14/h12-17H,1-11H2',
              'InChI=1S/C17H27Cl/c18-16-9-8-11-5-3-7-14-13-6-2-1-4-12(13)10-15(16)17(11)14/h11-17H,1-10H2',
              'InChI=1S/C4H10/c1-3-4-2/h3-4H2,1-2H3']
    print(Rollout.rollout_by_inchis(set([]), inchis, 10))


if __name__ == "__main__":
    test()
