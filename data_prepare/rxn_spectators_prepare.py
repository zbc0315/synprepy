# -*- coding: utf-8 -*-
# @Time     : 2021/4/25 21:48
# @Author   : zbc
# @Email    : zhangbc0315@outlook.com
# @FileName : rxn_spectators_prepare.py
# @Software : PyCharm

from typing import Iterator
import json

import pandas as pd
from rdkit.Chem import AllChem

from data_utils.known_spectators_loader import KnownSpectatorsLoader
from config import SSPath
from rmp_database.old_chem_db import OldChemDb


class RxnSpectatorsPrepare:

    _odb = OldChemDb()
    _inchi_to_smiles = {}

    _known_cats_smiles = []
    _known_sols_smiles = []

    @classmethod
    def _load_known_spectators(cls):
        cls._known_cats_smiles = list(KnownSpectatorsLoader.load_cat_smiles())
        cls._known_sols_smiles = list(KnownSpectatorsLoader.load_sol_smiles())

    @classmethod
    def _parse_sids(cls, sids: str) -> [int]:
        if sids is None or len(sids) == 0:
            return []
        else:
            return [int(sid) for sid in sids.split('.')]

    @classmethod
    def _get_all_rxn_before_2016(cls) -> Iterator:
        odb = OldChemDb()
        for rxn_data in odb.get_data_iter("public.clean_duplicate_rxn",
                                          ["rid", "rxn_code", "rxn_smi", "catalysts_sids", "solvents_sids", "tid_00"],
                                          None):
            rxn_data['catalysts_sids'] = cls._parse_sids(rxn_data['catalysts_sids'])
            rxn_data['solvents_sids'] = cls._parse_sids(rxn_data['solvents_sids'])
            yield rxn_data

    @classmethod
    def _smiles_from_inchi(cls, inchi: str) -> str:
        if inchi not in cls._inchi_to_smiles.keys():
            smiles = AllChem.MolToSmiles(AllChem.MolFromInchi(inchi))
            cls._inchi_to_smiles[inchi] = smiles
        return cls._inchi_to_smiles[inchi]

    @classmethod
    def _get_spectator_data_by_sid(cls, sid: int):
        spectator_data = cls._odb.search_one_data("public.spectators_mols",
                                                  ["inchi", "smiles", "been_catalyst_times", "been_solvent_times"],
                                                  f"sid={sid}")
        if spectator_data["smiles"] is None:
            spectator_data["smiles"] = cls._smiles_from_inchi(spectator_data["inchi"])
        return spectator_data

    @classmethod
    def _get_spectators_smiles_by_sids(cls, sids: [int], role: str):
        res = []
        for sid in sids:
            try:
                spectator_data = cls._get_spectator_data_by_sid(sid)
            except Exception as e:
                print(e)
                continue
            if role == "cat":
                times = spectator_data["been_catalyst_times"]
            else:
                times = spectator_data["been_solvent_times"]
            if times <= 1:
                continue
            res.append(spectator_data["smiles"])
        return res

    @classmethod
    def _is_mol_mapped(cls, mol) -> bool:
        for atom in mol.GetAtoms():
            if atom.GetAtomMapNum() is not None and atom.GetAtomMapNum() != 0:
                return True
        return False

    @classmethod
    def _clean_mols(cls, mols) -> ([str], [str]):
        """

        :param mols:
        :return: [mapped mol`s smiles cleaned map], [unmapped mol`s smiles]
        """
        mapped = []
        unmapped = []
        for mol in mols:
            if cls._is_mol_mapped(mol):
                mapped.append(AllChem.MolToSmiles(mol))
            else:
                unmapped.append(AllChem.MolToSmiles(mol))
        return mapped, unmapped

    @classmethod
    def _get_mol_role_by_known_spectators(cls, smiles: str):
        if smiles in cls._known_cats_smiles:
            return "cat"
        elif smiles in cls._known_sols_smiles:
            return "sol"
        return None

    @classmethod
    def _get_smiles_by_mid(cls, mid: int):
        mol_data = cls._odb.search_one_data("public.mols", ["inchi"], f"mid={mid}")
        return AllChem.MolToSmiles(AllChem.MolFromInchi(mol_data["inchi"]))

    @classmethod
    def _get_smiles_by_mids(cls, mids: [int]) -> [str]:
        return [cls._get_smiles_by_mid(mid) for mid in mids]

    @classmethod
    def _get_cats_sols_from_unmapped_reactants(cls, mapped_smiles: [str],
                                               mids: [int]) -> ([str], [str]):
        cats_smiles = []
        sols_smiles = []
        reactants_smiles = cls._get_smiles_by_mids(mids)
        for smiles in reactants_smiles:
            if smiles in mapped_smiles:
                continue
            role = cls._get_mol_role_by_known_spectators(smiles)
            if role == "cat":
                cats_smiles.append(smiles)
            elif role == "sol":
                sols_smiles.append(smiles)
        return cats_smiles, sols_smiles

    @classmethod
    def _clean_rxn_smiles(cls, rxn_smiles: str, rxn_code: str) -> (str, [str], [str]):
        rxn = AllChem.ReactionFromSmarts(rxn_smiles)
        products_smiles, unknown_smiles = cls._clean_mols(rxn.GetProducts())
        if len(unknown_smiles) != 0:
            raise ValueError("contain product with out mapped")
        reactants_smiles, unknown_smiles = cls._clean_mols(rxn.GetReactants())
        cats_smiles, sols_smiles = cls._get_cats_sols_from_unmapped_reactants(reactants_smiles,
                                                                              [int(mid) for mid in rxn_code.split(">")[0].split('.')])
        return reactants_smiles, products_smiles, cats_smiles, sols_smiles

    @classmethod
    def _extract_spectators(cls, rxn_data: {}):
        cats_smiles = cls._get_spectators_smiles_by_sids(rxn_data["catalysts_sids"], "cat")
        sols_smiles = cls._get_spectators_smiles_by_sids(rxn_data["solvents_sids"], "sol")
        reactants_smiles, products_smiles, cats_smiles_2, sols_smiles_2 = cls._clean_rxn_smiles(rxn_data["rxn_smi"],
                                                                                                rxn_data["rxn_code"])
        cats_smiles.extend(cats_smiles_2)
        sols_smiles.extend(sols_smiles_2)
        cats_smiles = list(set(cats_smiles))
        sols_smiles = list(set(sols_smiles))
        return reactants_smiles, products_smiles, cats_smiles, sols_smiles

    @classmethod
    def process_before_2016(cls):
        res_dic = {"rid": [], "tid": [], "reactants": [], "products": [], "cats": [], "sols": []}
        n = 0
        # e = 0
        for i, rxn_data in enumerate(cls._get_all_rxn_before_2016()):
            if i % 10000 == 0:
                # print(f"save {n}/{i} error: {e}")
                print(f"save {n}/{i}")
                df = pd.DataFrame(res_dic)
                df.to_csv(SSPath.RXN_WITH_SPECTATORS_FP, sep="\t", encoding="utf-8", index=False)
            rs_code = rxn_data["rxn_code"].split(">")[0]
            ps_code = rxn_data["rxn_code"].split(">")[-1]
            r_mids = rs_code.split('.')
            p_mids = ps_code.split('.')
            if any([p_mid in r_mids for p_mid in p_mids]):
                # e = e + 1
                continue
            cats = None
            sols = None
            try:
                reactants_smiles, products_smiles, cats_smiles, sols_smiles = cls._extract_spectators(rxn_data)
            except Exception as e:
                print("="*20)
                print(rxn_data["rid"])
                print(rxn_data["rxn_smi"])
                print(e)
                continue
            if len(cats_smiles) <= 1:
                cats = json.dumps(cats_smiles)
            if len(sols_smiles) <= 1:
                sols = json.dumps(sols_smiles)
            if cats is None and sols is None:
                continue
            n += 1
            res_dic["tid"].append(rxn_data["tid_00"])
            res_dic["rid"].append(rxn_data["rid"])
            res_dic["reactants"].append(json.dumps(reactants_smiles))
            res_dic["products"].append(json.dumps(products_smiles))
            res_dic["cats"].append(cats)
            res_dic["sols"].append(sols)
        df = pd.DataFrame(res_dic)
        df.to_csv(SSPath.RXN_WITH_SPECTATORS_FP, sep="\t", encoding="utf-8", index=False)

    @classmethod
    def process(cls):
        cls._load_known_spectators()
        cls.process_before_2016()


if __name__ == '__main__':
    RxnSpectatorsPrepare.process()
