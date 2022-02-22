# -*- coding: utf-8 -*-
# @Time     : 2021/6/15 15:20
# @Author   : zbc
# @Email    : zhangbc0315@outlook.com
# @FileName : prepare_wrong_rxn.py
# @Software : PyCharm

import sys
import os
import json
from datetime import datetime

from memory_profiler import profile
import pandas as pd
from rdkit.Chem import AllChem
from rdkit import rdBase
rdBase.DisableLog('rdApp.*')

# from config import RDPath
# from rmp_database.old_chem_db import OldChemDb
from chem_utils.rxn_mapper import RxnMapper

from crkitpy.molgraph.molecule import Molecule
from crkitpy.rxnhandler.rxn_template_apply import RxnTemplateApply
from crkitpy.inchiparse.inchi_parser import InchiParser
from crkitpy.inchiparse.inchi_writer import InchiWriter
from crkitpy.smilesparse.smiles_writer import SmilesWriter
from crkitpy.smatesparse.smates_parser import SmatesParser
from crkitpy.smilesparse.rdkit_smiles_parser import RdkitSmilesParser

# DATA_DP = "E:\\Data\\RxnMolPredictor\\RxnDeterminer"
DATA_DP = "C:\\Users\\zhang\\OneDrive\\Projects\\RXN-PREDICTOR\\Data\\RxnDeterminer"
NUM_PROCESS = 10
CLEANED_DUPLICATE_RXN_FP = os.path.join(DATA_DP, "rxndb_public_clean_duplicate_rxn.tsv")
RIGHT_RXN_FP = os.path.join(DATA_DP, "right_and_high_yield_rxn.tsv")
# WRONG_RXN_DP = "/home/hfnl/zbc/Projects/Data/wrong_rxn"
WRONG_RXN_FPS = os.path.join(DATA_DP, "wrong_rxn_{0}.tsv")
TEMP_FP = os.path.join(DATA_DP, "rxndb_public_rxn_template_00_filter.tsv")
# TEMP_FP = "/home/hfnl/zbc/Projects/Mcts/RxnTemplate/rxndb_public_rxn_template_00_filter.tsv"


class PrepareWrongRxn:

    _rxn_df = pd.read_csv(CLEANED_DUPLICATE_RXN_FP, sep='\t', encoding='utf-8')
    _right_and_high_yield_rxn_df = pd.read_csv(RIGHT_RXN_FP, sep='\t', encoding='utf-8')
    _temp_df = pd.read_csv(TEMP_FP, sep='\t', encoding='utf-8').query(f'rid_num>90')

    @classmethod
    def _get_temp(cls):
        """ 获得所有基础模板

        :return:
        """

        for _, row in cls._temp_df.sample(frac=1).iterrows():
            yield row['tid'], row['smates']
        # odb = OldChemDb()
        # for temp_data in odb.get_data_iter("public.rxn_template_00_filter", ["smates", "rid_num"], None, "rid_num", True):
        #     yield temp_data["smates"]

    @classmethod
    def _search_rxn_smi_from_rid(cls, rid: int):
        res = cls._rxn_df.query(f'rid=={rid}')
        return res.rxn_smi[res.index[0]]

    @classmethod
    def _get_right_rxn(cls):
        # odb = OldChemDb()
        for n, rid in enumerate(cls._right_and_high_yield_rxn_df['rid']):
            # rxn_data = odb.search_one_data('public.clean_duplicate_rxn', ['rxn_smi'], f'rid={rid}')
            rxn_smi = cls._search_rxn_smi_from_rid(rid)
            mapped_reactants, products_inchis = RxnMapper.get_mapped_reactants(rxn_smi)
            yield n, rid, mapped_reactants, products_inchis

    @classmethod
    def _rdmols_to_mols(cls, rdmols):
        return list(map(RdkitSmilesParser.rdmol_to_mol, rdmols))

    @classmethod
    def _is_saned_rxn(cls, products_inchis, reactants_inchis, right_product_inchi):
        """ 判断反应是否合理：
                产物中没有反应物则为合理

        :param products_inchis:
        :param reactants_inchis:
        :return:
        """
        if right_product_inchi in products_inchis:
            return False
        products_inchis_set = set(products_inchis)
        reactants_inchis_set = set(reactants_inchis)
        un_react_inchis = products_inchis_set.intersection(reactants_inchis_set)
        return len(un_react_inchis) == 0

    @classmethod
    def _get_main_product(cls, products: [Molecule]):
        main_mol: Molecule = products[0]
        if len(products) == 0:
            return main_mol

    @classmethod
    def _get_main_inchi(cls, inchis: [str]):
        inchis = sorted(inchis, key=lambda x: len(x))
        return inchis[-1]

    @classmethod
    def _get_wrong_rxn_by_temp(cls, right_rxn_reactants, smates, product_inchi, reactants_inchis):
        # res = []
        for products in RxnTemplateApply.reactants_to_products(smates, right_rxn_reactants, strict=True):
            try:
                ps_inchis = list(map(InchiWriter.mol_to_inchi, products))
            except Exception as e:
                # print(e)
                continue
            if not cls._is_saned_rxn(ps_inchis, reactants_inchis, product_inchi):
                continue
            # res.append(ps_inchis)
            yield ps_inchis
        # return res

    @classmethod
    def _get_wrong_rxn(cls, right_rxn_reactants_rd, product_inchi, rs_inchis):
        # res  = []
        for tid, temp_sma in cls._get_temp():
            right_rxn_reactants = cls._rdmols_to_mols(right_rxn_reactants_rd)
            for wrong_products_inchis in cls._get_wrong_rxn_by_temp(right_rxn_reactants, temp_sma, product_inchi, rs_inchis):
                yield cls._get_main_inchi(wrong_products_inchis), tid
                # res.append((cls._get_main_inchi(wrong_products_inchis), tid))
            # right_rxn_reactants.clear()
            # del right_rxn_reactants
            # del right_rxn_reactants
        # return res

    @classmethod
    def _write_wrong_rxns(cls, rid: int, rs_inchis: [str], p_inchi: str, wp_inchis: [str], out_fp: str, tids: [int]):
        if not os.path.exists(out_fp):
            with open(out_fp, 'w', encoding='utf-8')as f:
                f.write(f"rid\ttid\twp_inchi\trs_inchis\tp_inchi")
        with open(out_fp, 'a', encoding='utf-8')as f:
            for i, wp_inchi in enumerate(wp_inchis):
                f.write(f"\n{rid}\t{tids[i]}\t{wp_inchi}\t{json.dumps(rs_inchis)}\t{[p_inchi]}")

    @classmethod
    def _get_process_out_fp(cls, pid: int):
        return WRONG_RXN_FPS.format(pid)
        # return RDPath.WRONG_RXN_FPS.format(pid)

    @classmethod
    def _load_parsed_rids(cls, out_fp):
        if not os.path.exists(out_fp):
            return []
        df = pd.read_csv(out_fp, sep='\t', encoding='utf-8')
        return set(df['rid'])

    @classmethod
    def _is_mol_radical(cls, mol):
        if mol is None:
            return True
        for atom in mol.GetAtoms():
            if atom.GetNumRadicalElectrons() > 0:
                return True
        return False

    @classmethod
    def _show_tid_count(cls, tid_count):
        count = 0
        for key, value in tid_count.items():
            print(f'tid:{key} - num:{value}')
            count += value
        print(count)

    @classmethod
    def count_radical_mol(cls):
        res = {}
        for i in range(10):
            fp = cls._get_process_out_fp(i)
            df = pd.read_csv(fp, sep='\t', encoding='utf-8')
            for _, row in df.iterrows():
                mol = AllChem.MolFromInchi(row.wp_inchi)
                # if mol is None:
                #     print(row.wp_inchi)
                if not cls._is_mol_radical(mol):
                    continue
                tid = row.tid
                if tid not in res.keys():
                    cls._show_tid_count(res)
                    res[tid] = 0
                res[tid] += 1
            cls._show_tid_count(res)
        cls._show_tid_count(res)

    @classmethod
    def process(cls, pid):
        print(f'pid: {pid}')
        res_rxn = 16036
        ratio = 100
        needed_w_rxn = res_rxn * ratio
        out_fp = cls._get_process_out_fp(pid)
        parsed_rids = cls._load_parsed_rids(out_fp)
        wrong_idx = [177, 232, 392, 527, 1037, 1132, 1567, 1672, 1712,
                     1792, 1927, 2378, 9681, 5743, 1321, 733, 10231, 11165, 1,
                     185, 5669, 5679, 191, 420, 490, 439, 715, 1066, 834, 1101, 193]
        max_wrong_idx = [15070, 731, 14662, 543, 14664, 715, 586, 297, 428, 639]
        print(parsed_rids)
        for i, (lino, rid, reactants, products_inchis) in enumerate(cls._get_right_rxn()):
            if i <= max_wrong_idx[pid]:
                print(f"parsed rid: {i}-{rid}, max rid: {max_wrong_idx[pid]}")
                continue
            if rid in parsed_rids:
                print(f"parsed rid: {rid}")
                continue
            if lino % NUM_PROCESS != pid:
                continue
            # if i in wrong_idx:
            #     print(f"wrong i: {i} - {rid}")
            #     continue
            res_rxn -= 1
            print(f"{i} - rid:{rid} - {needed_w_rxn // res_rxn}")
            wps_inchis = []
            tids = []
            rs_inchis = list(map(AllChem.MolToInchi, reactants))
            if len(products_inchis) > 1:
                continue
            num = 0
            for n, (wrong_products_inchi, tid) in enumerate(cls._get_wrong_rxn(reactants, products_inchis[0], rs_inchis)):
                if wrong_products_inchi not in wps_inchis:
                    wps_inchis.append(wrong_products_inchi)
                    tids.append(tid)
                if n % ratio == 0:
                    print(f"wrong products inchis: {i} - {n} - {len(wps_inchis) + num} - {datetime.now()}")
                # if len(wps_inchis) >= needed_w_rxn // res_rxn:
                if False and len(wps_inchis) >= 1000:
                    num += len(wps_inchis)
                    cls._write_wrong_rxns(rid, rs_inchis, products_inchis[0], wps_inchis, out_fp, tids)
                    needed_w_rxn -= len(wps_inchis)
                    wps_inchis = []
                    tids = []
                    # break
                    # if num >= 10000:
                    #     break
            if len(wps_inchis) > 0:
                print(f'r{rid} has wrong rxn: {len(wps_inchis)}')
                cls._write_wrong_rxns(rid, rs_inchis, products_inchis[0], wps_inchis, out_fp, tids)
                needed_w_rxn -= len(wps_inchis)
            # break


if __name__ == '__main__':
    # print(sys.argv)
    # PrepareWrongRxn.process(4)
    PrepareWrongRxn.count_radical_mol()
