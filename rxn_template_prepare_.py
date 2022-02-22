#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Time    : 2022/1/5 23:56
# @Author  : zhangbc0315@outlook.com
# @File    : rxn_template_prepare_.py
# @Software: PyCharm

from tqdm import tqdm
import pandas as pd

from config import Config
from data_utils import RxnTemplates
from chem_utils.rxn_template.template_extractor import TemplateExtractor, RxnTemplateType


class RxnTemplatePrepare:

    # region ===== extract reaction templates =====

    @classmethod
    def _get_rxn_template(cls, rxn_smi: str, temp_type: RxnTemplateType):
        try:
            return TemplateExtractor.rxn_smiles_to_rxn_temp_smarts(rxn_smi, temp_type, True)
        except:
            return None

    @classmethod
    def _get_rows(cls, source_fp: str):
        rows = []
        with open(source_fp, 'r', encoding='utf-8')as f:
            for line in f.readlines():
                line = line.strip()
                if len(line) > 0:
                    rows.append(line)
        return rows

    @classmethod
    def _get_rid_rxn_smi_idx(cls, first_line: str) -> (int, int):
        columns = first_line.split('\t')
        if 'rxn_smi' not in columns:
            raise ValueError("Can't find column named 'rxn_smi'")
        if 'rid' not in columns:
            raise ValueError("Can't find column named 'rid")
        return columns.index('rid'), columns.index('rxn_smi')

    @classmethod
    def _init_result_file(cls, result_fp: str, first_line: str):
        with open(result_fp, 'w', encoding='utf-8')as f:
            f.write(f"{first_line}\n")

    @classmethod
    def _add_lines(cls, result_fp: str, lines: [str]):
        lines_block = '\n'.join(lines)
        with open(result_fp, 'a', encoding='utf-8')as f:
            f.write(f"{lines_block}\n")

    @classmethod
    def extract_rxn_templates(cls, config: Config):

        rid_idx = None
        rxn_smi_idx = None
        num_rct = 0
        num_ret = 0
        lines = cls._get_rows(config.rxn_data_tsv_file_path)
        results = []
        with tqdm(total=len(lines)) as pbar:
            for n, line in enumerate(lines):
                if n == 0:
                    rid_idx, rxn_smi_idx = cls._get_rid_rxn_smi_idx(line)
                    cls._init_result_file(config.rxn_code_with_rxn_template_tsv_file_path,
                                          "rid\trxn_centralized_template\trxn_extended_template")
                else:
                    rid = line.split('\t')[rid_idx]
                    rxn_smi = line.split('\t')[rxn_smi_idx]
                    rct = cls._get_rxn_template(rxn_smi, RxnTemplateType.CENTRALIZED)
                    ret = cls._get_rxn_template(rxn_smi, RxnTemplateType.EXTENDED)
                    if rct is not None:
                        num_rct += 1
                    if ret is not None:
                        num_ret += 1
                    results.append(f"{rid}\t{rct}\t{ret}")
                    if len(results) % 1000 == 0:
                        cls._add_lines(config.rxn_code_with_rxn_template_tsv_file_path, results)
                        results = []
                pbar.update(1)
                pbar.set_postfix_str(f"rct({num_rct}) ret({num_ret}) tot({n})")
        if len(results) > 0:
            cls._add_lines(config.rxn_code_with_rxn_template_tsv_file_path, results)

    # endregion

    # region ===== organize reaction template =====

    @classmethod
    def _organize_rxn_template(cls, templates_count: {str, int}, target_fp: str):
        templates_data = {'tid': [], 'rxn_template': [], 'num_of_covered_rxn': []}
        for i, (template, num) in enumerate(templates_count.items()):
            templates_data['tid'].append(i)
            templates_data['rxn_template'].append(template)
            templates_data['num_of_covered_rxn'].append(num)
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

    @classmethod
    def organize_rxn_template(cls, config: Config):
        rid_rxn_template_df = pd.read_csv(config.rxn_code_with_rxn_template_tsv_file_path, sep='\t', encoding='utf-8')
        sorted_centralized_templates_counts = cls._count_rxn_templates(rid_rxn_template_df.rxn_centralized_template)
        cls._organize_rxn_template(sorted_centralized_templates_counts, config.rxn_centralized_template_tsv_file_path)

        sorted_extended_templates_counts = cls._count_rxn_templates(rid_rxn_template_df.rxn_extended_template)
        cls._organize_rxn_template(sorted_extended_templates_counts, config.rxn_extended_template_tsv_file_path)

    # endregion

    # region ===== link rid and tid =====

    @classmethod
    def link_rid_and_tid(cls, config: Config):
        centralized_templates = RxnTemplates(config.rxn_centralized_template_tsv_file_path)
        extended_templates = RxnTemplates(config.rxn_extended_template_tsv_file_path)
        rid_with_temp_df = pd.read_csv(config.rxn_code_with_rxn_template_tsv_file_path, sep='\t', encoding='utf-8')
        idx = 0
        cls._init_result_file(config.rxn_code_with_rt_code_tsv_file_path, "rid\tc_tid\te_tid\tc_num\te_num")
        res_lines = []
        with tqdm(total=len(rid_with_temp_df), desc="link rid and tid")as pbar:
            for i, row in rid_with_temp_df.iterrows():
                c_temp = row.rxn_centralized_template
                e_temp = row.rxn_extended_template
                if pd.isna(c_temp) or c_temp == "None" or pd.isna(e_temp) or e_temp == "None":
                    continue
                rid = row.rid
                c_data = centralized_templates.query_by_temp(c_temp)
                e_data = extended_templates.query_by_temp(e_temp)
                res_lines.append(f"{rid}\t"
                                 f"{c_data.tid}\t"
                                 f"{e_data.tid}\t"
                                 f"{c_data.num_of_covered_rxn}\t"
                                 f"{e_data.num_of_covered_rxn}")
                if len(res_lines) == 10000:
                    cls._add_lines(config.rxn_code_with_rt_code_tsv_file_path, res_lines)
                    res_lines = []
                idx += 1
                pbar.update(1)
        if len(res_lines) > 0:
            cls._add_lines(config.rxn_code_with_rt_code_tsv_file_path, res_lines)
    # endregion

    @classmethod
    def process(cls):
        config = Config("config.json")
        cls.extract_rxn_templates(config)
        cls.organize_rxn_template(config)
        cls.link_rid_and_tid(config)


if __name__ == "__main__":
    RxnTemplatePrepare.process()
