#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Time    : 2021/8/13 10:50
# @Author  : zhangbc0315@outlook.com
# @File    : clean_rxndb.py
# @Software: PyCharm

from typing import Union
import os

from rdkit.Chem import AllChem
from rdkit import rdBase

from rmp_database.old_chem_db import OldChemDb
from config import SSPath

rdBase.DisableLog('rdApp.*')


class CleanRxndb:
    _source_tn = "public.clean_duplicate_rxn"
    _target_tn = "public.cleaned_rxn"
    _mols_tn = "public.mols"
    _cols = ['rid', 'rxn_code', 'rxn_smi', 'r_num', 'p_num', 'tid_10', 'tid_11', 'tid_20', 'tid_21', 'tid_00',
             'catalysts_sids', 'solvents_sids']
    _odb = OldChemDb()

    @classmethod
    def _get_rxn_data_from_source(cls):
        odb = OldChemDb()
        for rxn_data in odb.get_data_iter(cls._source_tn, cls._cols, None):
            yield rxn_data

    @classmethod
    def _add_rxn_data_to_target(cls, rxn_data):
        values = []
        cols = []
        for col in cls._cols:
            v = rxn_data[col]
            if v is None:
                continue
            values.append(v)
            cols.append(col)
        cls._odb.insert_one_data(cls._target_tn, cols, values, commit=False)

    @classmethod
    def _log_error_rxn(cls, rxn_data):
        if not os.path.exists(SSPath.CLEAN_RXN_ERROR_FP):
            with open(SSPath.CLEAN_RXN_ERROR_FP, 'w', encoding='utf-8') as f:
                f.write("\t".join(cls._cols))
        values = [str(rxn_data[col]) for col in cls._cols]
        with open(SSPath.CLEAN_RXN_ERROR_FP, 'a', encoding='utf-8')as f:
            f.write("\n")
            f.write("\t".join(values))

    # @classmethod
    # def _clean_mols_code(cls, mols_code: str) -> str:
    #     if len(mols_code) == 0:
    #         return mols_code
    #     mols_ids = mols_code.split('.')
    #     mols_ids = list(set(mols_ids))
    #     mols_ids = sorted(mols_ids, key=lambda x: int(x))
    #     return '.'.join(mols_ids)
    #
    # @classmethod
    # def clean_rxn_code(cls, rxn_code: str):
    #     rs_code, ss_code, ps_code = rxn_code.split('>')
    #     rs_code = cls._clean_mols_code(rs_code)
    #     ss_code = cls._clean_mols_code(ss_code)
    #     ps_code = cls._clean_mols_code(ps_code)
    #     return f"{rs_code}>{ss_code}>{ps_code}"

    @classmethod
    def _search_mol_id_by_inchi(cls, inchi: str) -> Union[int, None]:
        rxn_data = cls._odb.search_one_data(cls._mols_tn, ["mid"], f"inchi='{inchi}'")
        if rxn_data is not None:
            return rxn_data['mid']
        else:
            return None

    @classmethod
    def _get_mol_id_by_inchi(cls, inchi: str) -> int:
        mid = cls._search_mol_id_by_inchi(inchi)
        if mid is None:
            raise ValueError(f"inchi is not saved in mols: {inchi}")
        else:
            return mid

    @classmethod
    def _is_mapped(cls, mol):
        for atom in mol.GetAtoms():
            map_num = atom.GetAtomMapNum()
            if map_num is not None and map_num != 0:
                return True
        return False

    @classmethod
    def clean_rxn_smiles(cls, rxn_smi: str) -> (str, str):
        rxn = AllChem.ReactionFromSmarts(rxn_smi)
        res_rxn = AllChem.ChemicalReaction()
        res_rs_ids = set()
        res_ss_ids = set()
        res_ps_ids = set()
        for reactant in rxn.GetReactants():
            if cls._is_mapped(reactant):
                res_rxn.AddReactantTemplate(reactant)
                res_rs_ids.add(cls._get_mol_id_by_inchi(AllChem.MolToInchi(reactant)))
        for agent in rxn.GetAgents():
            if cls._is_mapped(agent):
                res_rxn.AddAgentTemplate(agent)
                res_ss_ids.add(cls._get_mol_id_by_inchi(AllChem.MolToInchi(agent)))
        for product in rxn.GetProducts():
            if cls._is_mapped(product):
                res_rxn.AddProductTemplate(product)
                res_ps_ids.add(cls._get_mol_id_by_inchi(AllChem.MolToInchi(product)))

        res_rs_ids = [str(i) for i in sorted(res_rs_ids, key=lambda x: x)]
        res_ss_ids = [str(i) for i in sorted(res_ss_ids, key=lambda x: x)] if len(res_ss_ids) > 0 else []
        res_ps_ids = [str(i) for i in sorted(res_ps_ids, key=lambda x: x)]

        return AllChem.ReactionToSmiles(
            res_rxn), f"{'.'.join(res_rs_ids)}>{'.'.join(res_ss_ids)}>{'.'.join(res_ps_ids)}"

    @classmethod
    def process(cls):
        for n, rxn_data in enumerate(cls._get_rxn_data_from_source()):
            try:
                rxn_smi, rxn_code = cls.clean_rxn_smiles(rxn_data["rxn_smi"])
            except:
                cls._log_error_rxn(rxn_data)
                continue
            rxn_data["rxn_smi"] = rxn_smi
            rxn_data["rxn_code"] = rxn_code
            cls._add_rxn_data_to_target(rxn_data)
            if n % 100 == 0 and n != 0:
                cls._odb.commit()
                print(n)
        cls._odb.commit()


if __name__ == "__main__":
    CleanRxndb.process()
