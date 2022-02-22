#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Time    : 2022/2/12 14:17
# @Author  : zhangbc0315@outlook.com
# @File    : rxn_template_prepare.py
# @Software: PyCharm

from tqdm import tqdm
import pandas as pd

from utils import TSVUtils
from config import Config
from data_utils import RxnData, RxnDataType, RxnTemplateData
from chem_utils.rxn_template.template_extractor import TemplateExtractor, RxnTemplateType


class RxnTemplatePrepare:

    def __init__(self):
        self._config = Config("config.json")
        self._rxn_data = RxnData(self._config.rxn_data_tsv_file_path,
                                 None,
                                 RxnDataType.TRAIN_TEST,
                                 self._config.evaluate_year)

    # region ===== extract rxn template =====

    @classmethod
    def _get_rxn_template(cls, rxn_smi: str, temp_type: RxnTemplateType):
        try:
            return TemplateExtractor.rxn_smiles_to_rxn_temp_smarts(rxn_smi, temp_type, True)
        except:
            return None

    def extract_rxn_templates(self):
        num_template = 0
        num_error = 0
        result_df = pd.DataFrame(columns=["rxn_code", "rct", "ret"])
        with tqdm(total=len(self._rxn_data)) as pbar:
            for n, rxn in enumerate(self._rxn_data.get_all_rxn(False)):
                rxn_code = rxn.rxn_code
                rxn_smiles = rxn.rxn_smiles
                rct = self._get_rxn_template(rxn_smiles, RxnTemplateType.CENTRALIZED)
                ret = self._get_rxn_template(rxn_smiles, RxnTemplateType.EXTENDED)
                if rct is None or ret is None:
                    num_error += 1
                    continue
                result_df = result_df.append({'rxn_code': rxn_code, 'rct': rct, 'ret': ret}, ignore_index=True)
                if len(result_df) >= 1000:
                    TSVUtils.df_to_tsv(result_df, self._config.rxn_code_with_rxn_template_tsv_file_path, mode='a')
                    result_df = pd.DataFrame(columns=["rxn_code", "rct", "ret"])
                num_template += 1
                pbar.update(1)
                pbar.set_postfix_str(f"template num: {num_template}, error num: {num_error}")
        if len(result_df) > 0:
            TSVUtils.df_to_tsv(result_df, self._config.rxn_code_with_rxn_template_tsv_file_path, mode='a')

    # endregion

    # region ===== filter rxn templates =====

    @classmethod
    def _organize_rxn_template(cls, templates_count: {str, int}, target_fp: str, pre: str):
        templates_data = {'rt_code': [], 'rxn_template': [], 'num_of_covered_rxn': [], 'rt_one_hot_idx': []}
        for i, (template, num) in enumerate(templates_count.items()):
            templates_data['rt_code'].append(f"{pre}{i}")
            templates_data['rxn_template'].append(template)
            templates_data['num_of_covered_rxn'].append(num)
            templates_data['rt_one_hot_idx'].append(i)
        df = pd.DataFrame(templates_data)
        df.to_csv(target_fp, sep='\t', encoding='utf-8', index=False)

    @classmethod
    def _count_rxn_templates(cls, all_templates: pd.Series):
        template_counts_series = all_templates.value_counts(sort=True, dropna=True)
        template_counts = {}
        for k, v in template_counts_series.items():
            if pd.isna(k) or k == "None" or k is None:
                continue
            template_counts[k] = v

        sorted_template_counts = {k: template_counts[k] for k in
                                  sorted(template_counts.keys(), key=lambda x: template_counts[x], reverse=True)}
        return sorted_template_counts

    def filter_rxn_templates(self):
        rxn_code_rt_df = pd.read_csv(self._config.rxn_code_with_rxn_template_tsv_file_path, sep='\t', encoding='utf-8')
        sorted_centralized_templates_counts = self._count_rxn_templates(rxn_code_rt_df.rct)
        self._organize_rxn_template(sorted_centralized_templates_counts, self._config.rxn_centralized_template_tsv_file_path, "TC")

        sorted_extended_templates_counts = self._count_rxn_templates(rxn_code_rt_df.ret)
        self._organize_rxn_template(sorted_extended_templates_counts, self._config.rxn_extended_template_tsv_file_path, "TE")

    # endregion

    def link_rxn_code_and_rt_code(self):
        rxn_code_rt_df = pd.read_csv(self._config.rxn_code_with_rxn_template_tsv_file_path, sep='\t', encoding='utf-8')
        rct_data = RxnTemplateData(self._config.rxn_centralized_template_tsv_file_path)
        ret_data = RxnTemplateData(self._config.rxn_extended_template_tsv_file_path)
        rxn_code_rt_code_df = pd.DataFrame(columns=["rxn_code", "rct_code", "ret_code", "c_num", "e_num"])
        with tqdm(total=len(rxn_code_rt_df))as pbar:
            pbar.set_description("link rxn code and rt code")
            for i, row in rxn_code_rt_df.iterrows():
                rxn_code = row.rxn_code
                rct_smiles = row.rct
                ret_smiles = row.ret
                rct = rct_data.get_template_by_rt_smiles(rct_smiles)
                ret = ret_data.get_template_by_rt_smiles(ret_smiles)
                rxn_code_rt_code_df = rxn_code_rt_code_df.append({"rxn_code": rxn_code,
                                                                  "rct_code": rct.rt_code,
                                                                  "ret_code": ret.rt_code,
                                                                  "c_num": rct.count,
                                                                  "e_num": ret.count}, ignore_index=True)
                if len(rxn_code_rt_code_df) >= 1000:
                    TSVUtils.df_to_tsv(rxn_code_rt_code_df, self._config.rxn_code_with_rt_code_tsv_file_path, 'a')
                    rxn_code_rt_code_df = pd.DataFrame(columns=["rxn_code", "rct_code", "ret_code", "c_num", "e_num"])
                pbar.update(1)
        if len(rxn_code_rt_code_df) > 0:
            TSVUtils.df_to_tsv(rxn_code_rt_code_df, self._config.rxn_code_with_rt_code_tsv_file_path, 'a')
        # rxn_code_rt_code_df.to_csv(self._config.rid_with_tid_tsv_file_path, sep='\t', encoding='utf-8', index=False)


if __name__ == "__main__":
    rtp = RxnTemplatePrepare()
    # rtp.extract_rxn_templates()
    rtp.filter_rxn_templates()
    rtp.link_rxn_code_and_rt_code()
