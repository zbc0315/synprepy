# -*- coding: utf-8 -*-
# @Time     : 2021/3/24 14:19
# @Author   : zbc
# @Email    : zhangbc0315@outlook.com
# @FileName : base_temp_analyzer.py
# @Software : PyCharm

from typing import Iterator, Tuple, Dict

import os

import json
import pandas as pd
from rdkit.Chem import AllChem
from rdkit import rdBase
rdBase.DisableLog('rdApp.error')


from database.old_chem_db import OldChemDb
from log_utils.logger import Logger
from config.path import BaseTempPath


class BaseTempAnalyzer:

    TOTAL_MOL = 1045932

    @classmethod
    def _get_product_inchis(cls) -> Iterator[str]:
        db = OldChemDb()
        for mol_data in db.get_data_iter('ml_data.train_temp_selector', ['mid', 'product_inchi'], None):
            yield mol_data['product_inchi']
        for mol_data in db.get_data_iter('ml_data.test_temp_selector', ['mid', 'product_inchi'], None):
            yield mol_data['product_inchi']

    @staticmethod
    def _get_mol_elements(inchi: str) -> Iterator[Tuple[int, str]]:
        mol = AllChem.MolFromInchi(inchi)
        if mol is None:
            Logger.info(f"error inchi: {inchi}")
            return None, None
        for atom in mol.GetAtoms():
            yield atom.GetAtomicNum(), atom.GetSymbol().capitalize(), atom.GetTotalValence(), atom.GetNumRadicalElectrons(), atom.GetFormalCharge()

    @classmethod
    def _get_all_elements(cls):
        atomic_nums = [1]
        elements = ['H']
        valences = [[1]]
        radicals = [[0]]
        charges = [[0]]
        for i, inchi in enumerate(cls._get_product_inchis()):
            if i % 10000 == 0 and i != 0:
                Logger.info(f"GetAllElements: process mol: {i}, get elements: {len(elements)}")
                # break
            for atomic_num, symbol, valence, radical, charge in cls._get_mol_elements(inchi):
                if atomic_num is None:
                    continue
                # atomic_num = str(atomic_num)
                if atomic_num not in atomic_nums:
                    atomic_nums.append(atomic_num)
                    elements.append(symbol)
                    valences.append([valence])
                    radicals.append([radical])
                    charges.append([charge])
                else:
                    idx = atomic_nums.index(atomic_num)
                    if valence not in valences[idx]:
                        valences[idx].append(valence)
                    if radical not in radicals[idx]:
                        radicals[idx].append(radical)
                    if charge not in charges[idx]:
                        charges[idx].append(charge)
        res = {}
        for i, num in enumerate(atomic_nums):
            res[num] = {'symbol': elements[i],
                        'valences': valences[i],
                        'radicals': radicals[i],
                        'charges': charges[i]}
        with open(BaseTempPath.ALL_ELEMENTS_FP, 'w', encoding='utf-8')as f:
            f.write(json.dumps(res))
        # res_atomic = ','.join(atomic_nums)
        # res_symbol = ','.join(elements)
        # res_valenc = ','.join(['.'.join(vs) for vs in valences])
        # res_radica = ','.join(['.'.join(rs) for rs in radicals])
        # res_charge = ','.join(['.'.join(cs) for cs in charges])
        # with open(BaseTempPath.ALL_ELEMENTS_FP, 'w', encoding='utf-8')as f:
        #     f.write(res_atomic + '\n' + res_symbol + '\n' + res_valenc + '\n' + res_radica + '\n' + res_charge)

    @classmethod
    def _process_mol(cls, inchi: str, elements: Dict):
        mol = AllChem.MolFromInchi(inchi)
        if mol is None:
            return None
        mol_analyzer = {}
        for key in elements.keys():
            atom_data = elements[key]
            symbol = atom_data['symbol']
            mol_analyzer[symbol] = 0
            for v in atom_data['valences']:
                mol_analyzer[f'{symbol}_V{v}'] = 0
            for r in atom_data['radicals']:
                mol_analyzer[f'{symbol}_R{r}'] = 0
            for c in atom_data['charges']:
                mol_analyzer[f'{symbol}_C{c}'] = 0
        mol_analyzer['Inchi'] = inchi
        mol_analyzer['Ring'] = 0
        mol_analyzer['Ring3'] = 0
        mol_analyzer['Ring4'] = 0
        mol_analyzer['Ring5'] = 0
        mol_analyzer['Ring6'] = 0
        mol_analyzer['Ring7'] = 0
        mol_analyzer['Ring8'] = 0
        mol_analyzer['Ring9~'] = 0
        mol_analyzer['RingAtom'] = 0
        mol_analyzer['Aromatic'] = 0
        mol_analyzer['AromaticAtom'] = 0
        mol_analyzer['Radical'] = 0

        for atom in mol.GetAtoms():

            mol_analyzer['H'] += atom.GetTotalNumHs()

            symbol = elements[str(atom.GetAtomicNum())]['symbol']
            mol_analyzer[symbol] += 1

            mol_analyzer[f'{symbol}_V{atom.GetTotalValence()}'] += 1
            mol_analyzer[f'{symbol}_R{atom.GetNumRadicalElectrons()}'] += 1
            mol_analyzer[f'{symbol}_C{atom.GetFormalCharge()}'] += 1

        ri = mol.GetRingInfo()
        mol_analyzer['Ring'] += ri.NumRings()
        for ring in ri.AtomRings():
            atom_num = len(ring)
            if atom_num < 9:
                mol_analyzer[f'Ring{atom_num}'] += 1
            else:
                mol_analyzer['Ring9~'] += 1
            mol_analyzer['RingAtom'] += atom_num
            if all([mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring]):
                mol_analyzer['Aromatic'] += 1
                mol_analyzer['AromaticAtom'] += atom_num
        return mol_analyzer

    @classmethod
    def _load_all_elements(cls):
        with open(BaseTempPath.ALL_ELEMENTS_FP, 'r', encoding='utf-8')as f:
            res = json.loads(f.read())
        return res

    @classmethod
    def _process_all_mol(cls, all_elements):
        mols_analyzer = None
        for n, inchi in enumerate(cls._get_product_inchis()):
            if n % 10000 == 0 and n != 0:
                print(n)
                # break
            mol_analyzer = cls._process_mol(inchi, all_elements)
            if mol_analyzer is None:
                continue
            if mols_analyzer is None:
                mols_analyzer = {}
                for key in mol_analyzer.keys():
                    mols_analyzer[key] = [mol_analyzer[key]]
            else:
                for key in mol_analyzer.keys():
                    mols_analyzer[key].append(mol_analyzer[key])
        df = pd.DataFrame(mols_analyzer)
        df.to_csv(BaseTempPath.ATOMS_COUNT_FP, sep='\t', encoding='utf-8', index=False)

    @classmethod
    def _count_contained_atom_mol_num(cls, element_count):
        res = 0
        for i in element_count.index:
            if i != 0:
                res += element_count[i]
        return res

    @classmethod
    def _analyze_atoms(cls, df, all_elements):
        elements = []
        for key in all_elements.keys():
            elements.append(all_elements[key]['symbol'])
        elements_res = {'element': [], 'contained_mol_num': [], 'contained_mol_ratio': []}
        for element in elements:
            Logger.info(f"analyze atoms: {element}")
            element_res = {'atom_num': [], 'mol_num': []}
            element_count = df.loc[:, element].value_counts()
            for i in range(element_count.index.min(), element_count.index.max() + 1):
                element_res['atom_num'].append(i)
                if i in element_count.index:
                    element_res['mol_num'].append(element_count[i])
                else:
                    element_res['mol_num'].append(0)
                element_df = pd.DataFrame(element_res)
                element_df.to_csv(os.path.join(BaseTempPath.RESULTS_DP, f'{element}.tsv'), sep='\t', encoding='utf-8', index=False)

            count = cls._count_contained_atom_mol_num(element_count)
            elements_res['element'].append(element)
            elements_res['contained_mol_num'].append(count)
            elements_res['contained_mol_ratio'].append(count/cls.TOTAL_MOL)
        elements_df = pd.DataFrame(elements_res)
        elements_df.to_csv(BaseTempPath.ELEMENTS_RESULT_FP, sep='\t', encoding='utf-8', index=False)

    @classmethod
    def _analyze_charge(cls, df, all_elements):
        charge_cols = []
        for key in all_elements.keys():
            symbol = all_elements[key]['symbol']
            for c in all_elements[key]['charges']:
                if c != 0:
                    charge_cols.append(f'{symbol}_C{c}')
        df['charge'] = 0
        for col in charge_cols:
            df['charge'] = df['charge'] + df[col]
        charge_count = df.loc[:, 'charge'].value_counts()
        num = cls.TOTAL_MOL - charge_count[0]
        ratio = num / cls.TOTAL_MOL
        with open(BaseTempPath.CHARGE_RESULT_FP, 'w', encoding='utf-8')as f:
            f.write(f'{num}({ratio}) mol has charge')

    @classmethod
    def _analyze_ring(cls, df, all_elements):
        cols = ['Ring', 'Ring9~']
        for n in range(3, 9):
            cols.append(f'Ring{n}')
        res = []
        for col in cols:
            ring_count_df = df.loc[:, col].value_counts()
            num = cls.TOTAL_MOL - ring_count_df[0]
            ratio = num / cls.TOTAL_MOL
            res.append(f'{num}({ratio}) mol has {col}')
        with open(BaseTempPath.RING_RESULT_FP, 'w', encoding='utf-8')as f:
            f.write('\n'.join(res))

    @classmethod
    def _analyze_aromatic(cls, df, all_elements):
        aromatic_count_df = df.loc[:, 'Aromatic'].value_counts()
        num = cls.TOTAL_MOL - aromatic_count_df[0]
        ratio = num / cls.TOTAL_MOL
        with open(BaseTempPath.AROMATIC_RESULT_FP, 'w', encoding='utf-8')as f:
            f.write(f'{num}({ratio}) mol has aromatic ring')

    @classmethod
    def _analyze_radical(cls, df, all_elements):
        cols = []
        for key in all_elements.keys():
            atom = all_elements[key]
            symbol = atom['symbol']
            for r in atom['radicals']:
                if r != 0:
                    cols.append(f'{symbol}_R{r}')
        df['Radical'] = 0
        for col in cols:
            df['Radical'] = df['Radical'] + df[col]
        radical_count_df = df.loc[:, 'Radical'].value_counts()
        num = cls.TOTAL_MOL - radical_count_df[0]
        ratio = num /cls.TOTAL_MOL
        with open(BaseTempPath.RADICAL_RESULT_FP, 'w', encoding='utf-8')as f:
            f.write(f'{num}({ratio}) mol has radical ring')

    @classmethod
    def _analyze_atoms_count(cls, all_elements):
        df = pd.read_csv(BaseTempPath.ATOMS_COUNT_FP, sep='\t', encoding='utf-8')
        # cls._analyze_atoms(df, all_elements)
        # cls._analyze_charge(df, all_elements)
        # cls._analyze_ring(df, all_elements)
        # cls._analyze_aromatic(df, all_elements)
        cls._analyze_radical(df, all_elements)

    @classmethod
    def process(cls):
        if not os.path.exists(BaseTempPath.ALL_ELEMENTS_FP):
            cls._get_all_elements()
        all_elements = cls._load_all_elements()
        if not os.path.exists(BaseTempPath.ATOMS_COUNT_FP):
            cls._process_all_mol(all_elements)
        cls._analyze_atoms_count(all_elements)


if __name__ == '__main__':
    BaseTempAnalyzer.process()
