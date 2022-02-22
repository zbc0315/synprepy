#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Time    : 2021/9/13 15:28
# @Author  : zhangbc0315@outlook.com
# @File    : mcts.py
# @Software: PyCharm
import os
from datetime import datetime

from mcts_lite.mol_utils import MolUtils
from mcts_lite.rxn_utils import RxnUtils
from mcts_lite.node_utils import NodeUtils
from mcts_lite.data.db_cac import DbCac
from mcts_lite.chem_utils import ChemUtils
from mcts_lite.policy.expand import Expand
from mcts_lite.policy.rollout import Rollout
from log_utils.logger import Logging
from mcts_lite.mcts_path import MctsPath


class Mcts:

    max_plan = 1
    max_deep = 25
    max_seconds = 300
    root_nid = 0
    current_nid = 0
    target_inchi = None

    @classmethod
    def init(cls, max_deep, max_seconds, max_plan):
        DbCac.init()

        cls.max_plan = max_plan
        cls.max_deep = max_deep
        cls.max_seconds = max_seconds
        cls.current_nid = cls.root_nid

    @classmethod
    def _init_target(cls, target_inchi):
        MolUtils.init()
        RxnUtils.init()
        NodeUtils.init(cls.max_deep)
        target_inchi = ChemUtils.remove_inchi_stereo(target_inchi)
        Logging.log.info(f'\n{"=" * 100}\nTarget inchi removed stereo: {target_inchi}')
        cls.target_inchi = target_inchi
        DbCac.target_inchi = target_inchi
        cls.root_nid = NodeUtils.create_root_node_by_inchi(target_inchi)

    @classmethod
    def _select_one_time(cls):
        cls.current_nid = cls.root_nid
        while not NodeUtils.is_leaf(cls.current_nid):
            child_nids_un_terminal = list(NodeUtils.get_un_terminal_child_nids(cls.current_nid))
            if len(child_nids_un_terminal) == 0:
                NodeUtils.set_no_select(cls.current_nid, True)
                return False
            child_nids_un_terminal_sorted = sorted(child_nids_un_terminal, key=lambda x: NodeUtils.get_A(x))
            cls.current_nid = child_nids_un_terminal_sorted[-1]
            NodeUtils.visit(cls.current_nid)
        return not NodeUtils.query_by_nid(cls.current_nid).no_expand

    @classmethod
    def select(cls):
        Logging.log.info('Select')
        while True:
            res = cls._select_one_time()
            if res:
                return True
            elif NodeUtils.query_by_nid(cls.root_nid).no_select:
                cls.current_nid = cls.root_nid
                return False

    @classmethod
    def expand(cls):
        Logging.log.info(f'Expand')
        Expand.expand_by_nid(cls.current_nid, cls.target_inchi)
        for child_nid in NodeUtils.query_by_nid(cls.current_nid).child_nids:
            child_node_data = NodeUtils.query_by_nid(child_nid)
            if child_node_data.is_solved:
                yield child_nid
        if cls.current_nid == cls.root_nid and NodeUtils.query_by_nid(cls.current_nid).no_select:
            Logging.log.warning(f'Can not expand from root: {cls.target_inchi}')

    @classmethod
    def rollout(cls):
        Logging.log.info('Rollout')
        for child_nid in NodeUtils.query_by_nid(cls.current_nid).child_nids:
            child_node_data = NodeUtils.query_by_nid(child_nid)
            if child_node_data.is_solved:
                reward = 2
            elif child_node_data.no_expand:
                reward = -1
            else:
                reward = Rollout.rollout_by_node(child_nid)
            NodeUtils.change_Q(reward/10, child_nid)
            NodeUtils.set_reward(child_nid, reward)

    @classmethod
    def update(cls):
        Logging.log.info('Update')
        parent_nid = cls.current_nid
        parent_node_data = NodeUtils.query_by_nid(parent_nid)
        while parent_node_data.level > 0:
            NodeUtils.update_Q(parent_nid)
            parent_nid = parent_node_data.parent_nid
            parent_node_data = NodeUtils.query_by_nid(parent_nid)

    @classmethod
    def search(cls, target_inchi):
        cls._init_target(target_inchi)
        begin_time = datetime.now()

        num_plan = 0
        search_num = 0

        while True:
            search_num += 1

            res = cls.select()
            if not res:
                yield cls.current_nid, False
                break

            for solved_nid in cls.expand():
                num_plan += 1
                yield solved_nid, True
                if num_plan >= cls.max_plan:
                    break
            if num_plan >= cls.max_plan:
                break

            used_seconds = (datetime.now() - begin_time).seconds
            Logging.log.info(f'used seconds: {used_seconds}')
            if used_seconds >= cls.max_seconds:
                cls.select()
                yield cls.current_nid, False
                break
            if search_num % 1 == 0:
                yield cls.current_nid, False

            cls.rollout()
            cls.update()


if __name__ == "__main__":
    Logging.init(999)
    Mcts.init(25, 30000, 1)
    # target = 'InChI=1S/C27H34FN3O5S/c1-26(2,3)36-25(32)31-15-13-30(14-16-31)22-18-20(29-37(33,34)23-8-5-4-7-21(23)28)17-19-9-12-27(10-6-11-27)35-24(19)22/h4-5,7-8,17-18,29H,6,9-16H2,1-3H3'
    # target = 'InChI=1S/C27H28N5O8P/c1-18-8-9-31(15-18)25-23-22(21(39-2)14-28-25)20(16-32(23)17-40-41(36,37)38)24(33)27(35)30-12-10-29(11-13-30)26(34)19-6-4-3-5-7-19/h3-9,14-16H,10-13,17H2,1-2H3,(H2,36,37,38)'
    # fes
    # target = 'InChI=1S/C25H26N7O8P/c1-16-27-14-32(28-16)23-21-20(19(39-2)12-26-23)18(13-31(21)15-40-41(36,37)38)22(33)25(35)30-10-8-29(9-11-30)24(34)17-6-4-3-5-7-17/h3-7,12-14H,8-11,15H2,1-2H3,(H2,36,37,38)'

    # ozanimod
    # target = 'InChI=1S/C23H24N4O3/c1-14(2)29-21-9-6-15(12-16(21)13-24)23-26-22(27-30-23)19-5-3-4-18-17(19)7-8-20(18)25-10-11-28/h3-6,9,12,14,20,25,28H,7-8,10-11H2,1-2H3'

    # remdesivir
    target = 'InChI=1S/C27H35N6O8P/c1-4-18(5-2)13-38-26(36)17(3)32-42(37,41-19-9-7-6-8-10-19)39-14-21-23(34)24(35)27(15-28,40-21)22-12-11-20-25(29)30-16-31-33(20)22/h6-12,16-18,21,23-24,34-35H,4-5,13-14H2,1-3H3,(H,32,37)(H2,29,30,31)/t17-,21+,23+,24+,27-,42-/m0/s1'

    # pf07321332
    # target = 'InChI=1S/C23H32F3N5O4/c1-21(2,3)16(30-20(35)23(24,25)26)19(34)31-10-13-14(22(13,4)5)15(31)18(33)29-12(9-27)8-11-6-7-28-17(11)32/h11-16H,6-8,10H2,1-5H3,(H,28,32)(H,29,33)(H,30,35)'

    # nifurtimox
    # target = 'InChI=1S/C10H13N3O5S/c1-8-7-19(16,17)5-4-12(8)11-6-9-2-3-10(18-9)13(14)15/h2-3,6,8H,4-5,7H2,1H3/b11-6+'
    # nifurtimox inter
    # target = 'InChI=1S/C5H12N2O2S/c1-5-4-10(8,9)3-2-7(5)6/h5H,2-4,6H2,1H3'

    # amisulpride
    target = 'InChI=1S/C17H27N3O4S/c1-4-20-8-6-7-12(20)11-19-17(21)13-9-16(25(22,23)5-2)14(18)10-15(13)24-3/h9-10,12H,4-8,11,18H2,1-3H3,(H,19,21)'

    # oliceridine
    target = 'InChI=1S/C22H30N2O2S/c1-25-18-7-15-27-19(18)16-23-13-10-21(20-6-2-5-12-24-20)11-14-26-22(17-21)8-3-4-9-22/h2,5-7,12,15,23H,3-4,8-11,13-14,16-17H2,1H3'

    # rimegepant
    target = 'InChI=1S/C28H28F2N6O3/c29-20-6-1-4-17(23(20)30)18-8-9-22(25-19(24(18)31)5-2-12-32-25)39-28(38)35-14-10-16(11-15-35)36-21-7-3-13-33-26(21)34-27(36)37/h1-7,12-13,16,18,22,24H,8-11,14-15,31H2,(H,33,34,37)'

    # ripretinib
    target = 'InChI=1S/C24H23BrFN5O2/c1-3-31-21-12-22(27-2)28-13-14(21)9-17(23(31)32)16-10-20(19(26)11-18(16)25)30-24(33)29-15-7-5-4-6-8-15/h4-8,10-13,17H,3,9H2,1-2H3,(H,27,28)(H2,29,30,33)'

    # risdiplam
    target = 'InChI=1S/C22H23N7O/c1-14-9-18(26-29-11-15(2)24-21(14)29)17-10-20(30)28-12-16(3-4-19(28)25-17)27-8-7-23-22(13-27)5-6-22/h3-4,9-12,23H,5-8,13H2,1-2H3'

    # tirbanibulin
    target = 'InChI=1S/C26H29N3O3/c30-26(28-19-21-4-2-1-3-5-21)18-24-9-6-23(20-27-24)22-7-10-25(11-8-22)32-17-14-29-12-15-31-16-13-29/h1-11,20H,12-19H2,(H,28,30)'

    # Berotralstat
    target = 'InChI=1S/C30H26F4N6O/c31-24-10-9-22(28(37-17-18-7-8-18)21-5-1-3-19(11-21)15-35)13-25(24)38-29(41)26-14-27(30(32,33)34)39-40(26)23-6-2-4-20(12-23)16-36/h1-6,9-14,18,28,37H,7-8,16-17,36H2,(H,38,41)'

    # Avapritinib
    # target = 'InChI=1S/C26H27FN10/c1-26(28,20-3-5-22(27)6-4-20)21-13-29-25(30-14-21)36-9-7-35(8-10-36)24-23-11-18(16-37(23)33-17-31-24)19-12-32-34(2)15-19/h3-6,11-17H,7-10,28H2,1-2H3/t26-/m0/s1'

    # Avapritinib inter
    # target = 'InChI=1S/C25H22FN9O/c1-32-14-20(13-30-32)18-10-22-24(29-16-31-35(22)15-18)33-6-8-34(9-7-33)25-27-11-19(12-28-25)23(36)17-2-4-21(26)5-3-17/h2-5,10-16H,6-9H2,1H3'
    for i, is_right in Mcts.search(target):
        if not is_right:
            continue
        if i != 0:
            print("="*50)
            RxnUtils.save(os.path.join(MctsPath.DEBUG_DP, f'{i}-Rxn.tsv'))
            NodeUtils.save(os.path.join(MctsPath.DEBUG_DP, f'{i}-Node.tsv'))
            MolUtils.save(os.path.join(MctsPath.DEBUG_DP, f'{i}-Mol.tsv'))
            smileses, _, __ = NodeUtils.get_rxn_plan(i)
            for smiles in smileses:
                print(smiles)
            # print(NodeUtils.get_rxn_plan(i))
