#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Time    : 2021/8/31 15:13
# @Author  : zhangbc0315@outlook.com
# @File    : mcts_path.py
# @Software: PyCharm

import os
import platform


class MctsPath:

    if platform.system() == 'Windows':
        _ROOT = 'C:\\Users\\zhang\\Documents\\Data\\RxnMolPredictor\\Mcts'
    else:
        _ROOT = '/home/hfnl/zbc/Projects/Mcts'
    # _ROOT = '/Users/zbc/SynologyDrive/Projects/RXN-PREDICTOR/Data/Mcts'
    LOG_DP = os.path.join(_ROOT, 'logs')
    LOG_FP = os.path.join(LOG_DP, 'log-%d.txt')
    DEBUG_DP = os.path.join(_ROOT, 'debug')
    BASIC_TIDS_FP = os.path.join(_ROOT, 'tids')
    COMP_TIDS_FP = os.path.join(_ROOT, 'comp_tids')
    SPEC_TIDS_FP = os.path.join(_ROOT, 'spectators_tids.txt')

    # =======================================
    # -------------     CAC     -------------

    _CAC_DP = os.path.join(_ROOT, 'CAC')

    CAC_INCHI_BBM_FP = os.path.join(_CAC_DP, 'rxndb_public_building_block_mols.tsv')
    CAC_INCHI_SM_FP = os.path.join(_CAC_DP, 'rxndb_public_for_sale_mols.tsv')
    CAC_INCHI_RM_FP = os.path.join(_CAC_DP, 'rxndb_public_reactant_mols.tsv')

    CAC_INCHI_NON_STEREO_BBM_FP = os.path.join(_CAC_DP, 'non_stereo_rxndb_public_building_block_mols.tsv')
    CAC_INCHI_NON_STEREO_SM_FP = os.path.join(_CAC_DP, 'non_stereo_rxndb_public_for_sale_mols.tsv')
    CAC_INCHI_NON_STEREO_RM_FP = os.path.join(_CAC_DP, 'non_stereo_rxndb_public_reactant_mols.tsv')

    CAC_SMILES_BBM_FP = os.path.join(_CAC_DP, 'smiles_rxndb_public_building_block_mols.tsv')
    CAC_SMILES_SM_FP = os.path.join(_CAC_DP, 'smiles_rxndb_public_for_sale_mols.tsv')
    CAC_SMILES_RM_FP = os.path.join(_CAC_DP, 'smiles_rxndb_public_reactant_mols.tsv')

    CAC_SMILES_NON_STEREO_BBM_FP = os.path.join(_CAC_DP, 'smiles_non_stereo_rxndb_public_building_block_mols.tsv')
    CAC_SMILES_NON_STEREO_SM_FP = os.path.join(_CAC_DP, 'smiles_non_stereo_rxndb_public_for_sale_mols.tsv')
    CAC_SMILES_NON_STEREO_RM_FP = os.path.join(_CAC_DP, 'smiles_non_stereo_rxndb_public_reactant_mols.tsv')

    # =====================

    _MODEL_DP = os.path.join(_ROOT, 'Models')

    # ============================
    # --  Model Basic Template  --

    MODEL_BT_FP = os.path.join(_MODEL_DP, 'bt_model_200_2.185096339638548_0.39687987110802564.pt')

    # ============================
    # -- Model Complex Template --

    MODEL_CT_FP = os.path.join(_MODEL_DP, 'ct_model_141_2.0868713010085855_0.457849478390462_.pt')

    # ===============================
    # -- Model Rxn Sane Determiner --

    MODEL_RD_FP = os.path.join(_MODEL_DP, 'rd_model_opt_74_a0.9702884915741274_p0.9720672494896699_r0.9409673341646233_.pt')
    MODEL_RD_WITHOUT_CS_FP = os.path.join(_MODEL_DP, 'rd_without_cat_sol_model_opt_51.pt')

    # ===========================
    # --     Data Template     --

    _TEMP_DP = os.path.join(_ROOT, 'RxnTemplate')

    BASIC_RXN_TEMP_FP = os.path.join(_TEMP_DP, 'rxndb_public_rxn_template_00_filter.tsv')
    COMP_RXN_TEMP_FP = os.path.join(_TEMP_DP, 'rxndb_public_rxn_template_10_filter.tsv')

    # ==========================
    # --       Catalyst       --

    _CAT_DP = os.path.join(_ROOT, 'CatPredictor')
    CATS_COUNT_FILTER_FP = os.path.join(_CAT_DP, 'cats_count_filter.tsv')
    CATS_MODEL_FP = os.path.join(_CAT_DP, 'model_192_l0.004000192781760229_c0.5485781990521327.pt')
    CATS_IDX_TO_WEBMID_FP = os.path.join(_CAT_DP, 'cats_idx_to_webmid.tsv')

    # =========================
    # --       Solvent       --

    _SOL_DP = os.path.join(_ROOT, 'SolPredictor')
    SOLS_COUNT_FILTER_FP = os.path.join(_SOL_DP, 'sols_count_filter.tsv')
    SOLS_MODEL_FP = os.path.join(_SOL_DP, 'model_98_l0.004448515211598848_c0.4039088460872167.pt')
    SOLS_IDX_TO_WEBMID_FP = os.path.join(_SOL_DP, 'sols_idx_to_webmid.tsv')

    # =============================
    # --      Mcts Evaluate      --

    _EVA_DP = os.path.join(_ROOT, 'Evaluate')
    EVA10_DP = os.path.join(_ROOT, 'Evaluate10')
    EVA30_DP = os.path.join(_ROOT, 'Evaluate30')
    TARGET_PRODUCTS_FP = os.path.join(_EVA_DP, 'target_products.tsv')
    EVALUATE_RESULT_FP = os.path.join(_EVA_DP, 'evaluate_result_%d.tsv')
    EVALUATE30_RESULT_FP = os.path.join(EVA30_DP, 'evaluate_result_%d.tsv')
    ALL_EXAMPLE_DP = os.path.join(_ROOT, 'all_examples')
    ALL_EXAMPLE_FP = os.path.join(_ROOT, 'all_examples.docx')


if __name__ == "__main__":
    pass
