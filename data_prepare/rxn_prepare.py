# -*- coding: utf-8 -*-
# @Time     : 2021/3/30 16:25
# @Author   : zbc
# @Email    : zhangbc0315@outlook.com
# @FileName : rxn_prepare.py
# @Software : PyCharm

from typing import Dict, Union, Set
import json

from rdkit.Chem import AllChem

from log_utils.logger import Logger
from rmp_database.old_chem_db import OldChemDb
from rmp_database.chem_db_orm import Session, MolTable, RxnTable, RxnMolTable


class RxnPrepare:

    """ 从原始数据中，重新清洗数据，存储到新数据库

    """

    _ocdb = OldChemDb()
    _saved_mols_smiles = set()
    _session = Session()
    _unsaved_mols: [MolTable] = []

    # ================
    # -- mol tables --

    @classmethod
    def _get_saved_mol_smiles(cls) -> None:
        """ 获得已经存储并且已经commit的所有分子的smiles

        :return:
        """
        for mol in cls._session.query(MolTable):
            cls._saved_mols_smiles.add(mol.smiles)
        Logger.info(f"Get Saved Mol Smiles: {len(cls._saved_mols_smiles)}")

    @classmethod
    def _get_uncommitted_mol_by_smiles(cls, smiles: str) -> MolTable:
        """ 根据smiles，获得已经存储，但是未commit的分子

        :param smiles:
        :return:
        """
        for mol in cls._unsaved_mols:
            if mol.smiles == smiles:
                return mol

    @classmethod
    def _get_saved_mol_by_smiles(cls, smiles: str) -> MolTable:
        """ 根据smiles，获得已经存储并且已经commit的分子

        :param smiles:
        :return:
        """
        mol = cls._session.query(MolTable).filter(MolTable.smiles == smiles).one()
        return mol

    # =======================
    # -- process 2016 2018 --

    @classmethod
    def process_16_18(cls):
        pass

    # =======================
    # -- process 2001 2016 --

    @classmethod
    def _get_origin_rxn_by_rid(cls, rid: int) -> Union[Dict, None]:
        """ 通过rid，从老数据库中获得化学反应的来源信息

        :param rid:
        :return:
        """
        origin_rxn_data = cls._ocdb.search_one_data('public.origin_rxn', ['source'], f'id={rid}')
        return origin_rxn_data

    @classmethod
    def _get_spectator_inchi_by_sid(cls, sid: int) -> str:
        """ 通过sid，从老数据库中获得催化剂、溶剂的Inchi

        :param sid:
        :return:
        """
        spectator_data = cls._ocdb.search_one_data('public.spectators_mols', ['inchi'], f'sid={sid}')
        return spectator_data['inchi']

    @classmethod
    def _parse_spectators_sids(cls, sids: str) -> [str]:
        """ 解析形如'454.4165'的字符串到Inchi的列表，如果sids为None，则传出空列表

        :param sids:
        :return:
        """
        if sids is None or len(sids) == 0:
            return []
        else:
            return [cls._get_spectator_inchi_by_sid(int(sid)) for sid in sids.split('.')]

    @classmethod
    def _clear_map_num(cls, mol):
        """ 清除分子的atom map num

        :param mol:
        :return:
        """
        for atom in mol.GetAtoms():
            atom.SetAtomMapNum(0)
        return mol

    @classmethod
    def _parse_rxn_smi(cls, rxn_smi: str) -> ([str], [str]):
        """ 解析反应Smiles，获得反应物的Smiles列表，产物的Smiles列表

        :param rxn_smi:
        :return:
        """
        rxn = AllChem.ReactionFromSmarts(rxn_smi)
        rs = []
        ps = []
        for r in rxn.GetReactants():
            r = cls._clear_map_num(r)
            rs.append(AllChem.MolToSmiles(r))
        for p in rxn.GetProducts():
            p = cls._clear_map_num(p)
            ps.append(AllChem.MolToSmiles(p))
        return rs, ps

    @classmethod
    def _parse_source(cls, source: str) -> (str, str, str, int):
        """ 解析反应来源

        :param source:
        :return: 文档id，标题，段落，年份
        """
        source_json = json.loads(source)
        return source_json['doc_id'], source_json['h_text'], source_json['p_text'], int(source_json['doc_id'][2:6])

    @classmethod
    def _save_mols_by_rdmol(cls, mol) -> MolTable:
        """ 存储mol到数据库

        :param mol:
        :return:
        """
        smiles = AllChem.MolToSmiles(mol)
        if smiles in cls._saved_mols_smiles:
            mol = cls._get_uncommitted_mol_by_smiles(smiles)
            if mol is None:
                mol = cls._get_saved_mol_by_smiles(smiles)
            return mol
        else:
            inchi = AllChem.MolToInchi(mol)
            mol = MolTable(mid=len(cls._saved_mols_smiles),
                           smiles=smiles,
                           inchi=inchi)
            cls._session.add(mol)
            cls._saved_mols_smiles.add(smiles)
            cls._unsaved_mols.append(mol)
            return mol

    @classmethod
    def _save_mols_by_smiles(cls, smileses: [str]) -> [MolTable]:
        return [cls._save_mols_by_rdmol(AllChem.MolFromSmiles(smiles)) for smiles in smileses]

    @classmethod
    def _save_mols_by_inchis(cls, inchises: [str]) -> [MolTable]:
        return [cls._save_mols_by_rdmol(AllChem.MolFromInchi(inchi)) for inchi in inchises]

    @classmethod
    def _get_rxn_code(cls, rs: [MolTable], ps: [MolTable]) -> str:
        rids = '.'.join([str(r.mid) for r in rs])
        pids = '.'.join([str(p.mid) for p in ps])
        return f'{rids}>>{pids}'

    @classmethod
    def _save_rxn_mol_relations_by_role(cls, rxn: RxnTable, mols: [MolTable], role: int) -> [RxnMolTable]:
        """

        :param rxn:
        :param mols:
        :param role: 0-反应物；2-催化剂；3-溶剂；1-产物
        :return:
        """
        res = []
        for mol in mols:
            rxn_mol = RxnMolTable(mol_role=role, rid=rxn.rid, mid=mol.mid)
            cls._session.add(rxn_mol)
            # mol_rxns = mol.rxn_mols
            # mol_rxns.append(rxn_mol)
            # mol.rxn_mols = mol_rxns
            res.append(rxn_mol)
        return res

    @classmethod
    def _save_rxn_mol_relations(cls, rxn: RxnTable, rs: [MolTable], ps: [MolTable],
                                cats: [MolTable], sols: [MolTable]):
        rxn_rs = cls._save_rxn_mol_relations_by_role(rxn, rs, 0)
        rxn_ps = cls._save_rxn_mol_relations_by_role(rxn, ps, 1)
        rxn_cats = cls._save_rxn_mol_relations_by_role(rxn, cats, 2)
        rxn_sols = cls._save_rxn_mol_relations_by_role(rxn, sols, 3)
        rxn.rxn_mols = rxn_rs + rxn_ps + rxn_cats + rxn_sols

    @classmethod
    def process_01_16(cls):
        ocdb = OldChemDb()
        for rxn_data in ocdb.get_data_iter('public.clean_duplicate_rxn',
                                           ['rid', 'rxn_code', 'rxn_smi', 'r_num', 'p_num', 'catalysts_sids', 'solvents_sids'],
                                           None):
            origin_rxn_data = cls._get_origin_rxn_by_rid(rxn_data['rid'])
            rid = rxn_data['rid']
            rxn_smi = rxn_data['rxn_smi']
            r_smileses, p_smileses = cls._parse_rxn_smi(rxn_smi)
            cat_inchises = cls._parse_spectators_sids(rxn_data['catalysts_sids'])
            sol_inchises = cls._parse_spectators_sids(rxn_data['solvents_sids'])
            rs = cls._save_mols_by_smiles(r_smileses)
            ps = cls._save_mols_by_smiles(p_smileses)
            cats = cls._save_mols_by_inchis(cat_inchises)
            sols = cls._save_mols_by_inchis(sol_inchises)

            source, title, para, year = cls._parse_source(origin_rxn_data['source'])
            rxn_code = cls._get_rxn_code(rs, ps)
            rxn = RxnTable(rid=rid, code=rxn_code, rxn_smiles=rxn_smi, source=source,
                           title=title, paragraph=para, year=year)
            cls._save_rxn_mol_relations(rxn, rs, ps, cats, sols)
            cls._session.add(rxn)
            cls._session.commit()

    @classmethod
    def process(cls):
        cls._get_saved_mol_smiles()
        cls.process_01_16()
        cls.process_16_18()


if __name__ == '__main__':
    RxnPrepare.process()
