# -*- coding: utf-8 -*-
# @Time     : 2021/3/24 14:20
# @Author   : zbc
# @Email    : zhangbc0315@outlook.com
# @FileName : mcts_path.py
# @Software : PyCharm

import os

# _PROJ_PATH = "C:/Users/zhang/Documents/Data/RxnMolPredictor"
_PROJ_PATH = "E:\\Data\\RxnMolPredictor"


class Path:

    dirs = []

    @classmethod
    def init(cls):
        return
        for dp in cls.dirs:
            if not os.path.exists(dp):
                os.mkdir(dp)


class BaseTempPath(Path):

    _PATH = os.path.join(_PROJ_PATH, "BaseTemp")
    _DATA_ANALYZER_DP = os.path.join(_PATH, "DataAnalyzer")

    ALL_ELEMENTS_FP = os.path.join(_DATA_ANALYZER_DP, 'all_elements.json')
    ATOMS_COUNT_FP = os.path.join(_DATA_ANALYZER_DP, 'atoms_count.tsv')

    RESULTS_DP = os.path.join(_DATA_ANALYZER_DP, 'results')
    ELEMENTS_RESULT_FP = os.path.join(RESULTS_DP, 'elements.tsv')
    CHARGE_RESULT_FP = os.path.join(RESULTS_DP, 'charge.txt')
    RING_RESULT_FP = os.path.join(RESULTS_DP, 'ring.txt')
    AROMATIC_RESULT_FP = os.path.join(RESULTS_DP, 'aromatic.txt')
    RADICAL_RESULT_FP = os.path.join(RESULTS_DP, 'radical.txt')

    _HIGHWAY_GNN_DP = os.path.join(_PATH, 'highway_gnn')
    TRAIN_BATCH_DATA_DP = os.path.join(_HIGHWAY_GNN_DP, 'train_batch_data_1000')
    TEST_BATCH_DATA_DP = os.path.join(_HIGHWAY_GNN_DP, 'test_batch_data_1000')
    HG_MODELS_DP = os.path.join(_HIGHWAY_GNN_DP, 'models')
    TIDS_FP = os.path.join(_PATH, "tids")

    CBE_DP = os.path.join(_PATH, "consider_bond_energy")
    CBE_SINGLE_DATA_TRAIN_DP = os.path.join(CBE_DP, "train_single_data")
    CBE_SINGLE_DATA_TEST_DP = os.path.join(CBE_DP, "test_single_data")

    CBE_BATCH_DATA_TRAIN_DP = os.path.join(CBE_DP, "train_batch_1000_data")
    CBE_BATCH_DATA_TEST_DP = os.path.join(CBE_DP, "test_batch_1000_data")

    CBE_MODEL_DP = os.path.join(CBE_DP, "model")
    CBE_TIDS_FP = os.path.join(CBE_DP, "cbe_tids.txt")
    CBE_NON_RD_SMILES_FP = os.path.join(CBE_DP, "non_rdkit_smiles.txt")
    CBE_NON_BDE_SMILES_FP = os.path.join(CBE_DP, "non_bde_smiles.txt")

    BOND_ENERGY_MODEL = os.path.join(_PATH, "bond_energy.pt")

    dirs = [_PATH, _DATA_ANALYZER_DP, RESULTS_DP, _HIGHWAY_GNN_DP, TRAIN_BATCH_DATA_DP, TEST_BATCH_DATA_DP,
            HG_MODELS_DP, CBE_DP, CBE_SINGLE_DATA_TRAIN_DP, CBE_SINGLE_DATA_TEST_DP,
            CBE_BATCH_DATA_TRAIN_DP, CBE_BATCH_DATA_TEST_DP, CBE_MODEL_DP]


class CompTempPath(Path):

    _PATH = os.path.join(_PROJ_PATH, "CompTemp")

    COMP_TIDS_FP = os.path.join(_PATH, 'comp_tids')

    _GCN_DP = os.path.join(_PATH, 'gcn')
    GCN_SINGLE_DATA_DP = os.path.join(_GCN_DP, 'single_data')
    GCN_TRAIN_BATCH_DATA_DP = os.path.join(_GCN_DP, 'train_batch_1000_data')
    GCN_TEST_BATCH_DATA_DP = os.path.join(_GCN_DP, 'test_batch_1000_data')
    MODEL_DP = os.path.join(_PATH, 'models')

    _HIGHWAY_GCN_DP = os.path.join(_PATH, 'highway_gcn')
    HIGHWAY_TRAIN_BATCH_DATA_DP = os.path.join(_HIGHWAY_GCN_DP, "train_batch_1000_data")
    HIGHWAY_TEST_BATCH_DATA_DP = os.path.join(_HIGHWAY_GCN_DP, "test_batch_1000_data")

    GCN_RESULT_FP = os.path.join(_GCN_DP, "gcn_result.tsv")
    HIGHWAY_RESULT_FP = os.path.join(_HIGHWAY_GCN_DP, "highway_result.tsv")

    dirs = [_PATH, _GCN_DP, GCN_SINGLE_DATA_DP, _HIGHWAY_GCN_DP, GCN_TRAIN_BATCH_DATA_DP, GCN_TEST_BATCH_DATA_DP, MODEL_DP]


class SpectatorsSelectorPath(Path):

    # =========================
    # -- Spectators Selector --

    _SS_APP_DP = os.path.join(_PROJ_PATH, "Spectators Selector")
    SS_NOTE_FP = os.path.join(_PROJ_PATH, "note.txt")
    ALL_CATS_SIDS_FP = os.path.join(_SS_APP_DP, "all_cats_sids.txt")
    ALL_SOLS_SIDS_FP = os.path.join(_SS_APP_DP, "all_sols_sids.txt")

    KNOWN_CAT_INCHIS = os.path.join(_SS_APP_DP, "Spectators/knownCatalystInChIs.txt")
    KNOWN_SOL_INCHIS = os.path.join(_SS_APP_DP, "Spectators/knownSolventInChIs.txt")

    RXN_WITH_SPECTATORS_FP = os.path.join(_SS_APP_DP, "rxn_with_spectators_before_2016.tsv")
    TIDS_FP = os.path.join(_SS_APP_DP, "spectators_tids.txt")

    CATS_COUNT_FP = os.path.join(_SS_APP_DP, "cats_count.tsv")
    SOLS_COUNT_FP = os.path.join(_SS_APP_DP, "sols_count.tsv")
    CATS_COUNT_WITH_SID_FP = os.path.join(_SS_APP_DP, "cats_count_with_sid.tsv")
    SOLS_COUNT_WITH_SID_FP = os.path.join(_SS_APP_DP, "sols_count_with_sid.tsv")

    CATS_TIMES_COUNT_FP = os.path.join(_SS_APP_DP, "cats_times_count.tsv")
    SOLS_TIMES_COUNT_FP = os.path.join(_SS_APP_DP, "sols_times_count.tsv")

    CATS_COUNT_FILTER_FP = os.path.join(_SS_APP_DP, "cats_count_filter.tsv")
    SOLS_COUNT_FILTER_FP = os.path.join(_SS_APP_DP, "sols_count_filter.tsv")

    _CAT_GNN_DATA_DP = os.path.join(_SS_APP_DP, "cat_gnn_data")
    _SOL_GNN_DATA_DP = os.path.join(_SS_APP_DP, "sol_gnn_data")

    _CAT_GNN_DATA_DP = "C:\\Users\\zhang\\Documents\\Data\\RxnMolPredictor\\SpectatorsSelector\\cat"
    _SOL_GNN_DATA_DP = "C:\\Users\\zhang\\Documents\\Data\\RxnMolPredictor\\SpectatorsSelector\\sol"

    CAT_GNN_SINGLE_DATA_DP = os.path.join(_CAT_GNN_DATA_DP, "gnn_single_data")
    SOL_GNN_SINGLE_DATA_DP = os.path.join(_SOL_GNN_DATA_DP, "gnn_single_data")

    TRAIN_CAT_GNN_BATCH_DATA_DP = os.path.join(_CAT_GNN_DATA_DP, "train_gnn_batch_500_data")
    TRAIN_SOL_GNN_BATCH_DATA_DP = os.path.join(_SOL_GNN_DATA_DP, "train_gnn_batch_500_data")

    TEST_CAT_GNN_BATCH_DATA_DP = os.path.join(_CAT_GNN_DATA_DP, "test_gnn_batch_500_data")
    TEST_SOL_GNN_BATCH_DATA_DP = os.path.join(_SOL_GNN_DATA_DP, "test_gnn_batch_500_data")

    CAT_GNN_DATA_INFO_FP = os.path.join(_CAT_GNN_DATA_DP, "cat_gnn_single_data_info.tsv")
    SOL_GNN_DATA_INFO_FP = os.path.join(_SOL_GNN_DATA_DP, "sol_gnn_single_data_info.tsv")
    CAT_TID_COR_FP = os.path.join(_CAT_GNN_DATA_DP, "cat_tid_cor.tsv")
    SOL_TID_COR_FP = os.path.join(_SOL_GNN_DATA_DP, "sol_tid_cor.tsv")

    CAT_GNN_TOPK_FP = os.path.join(_CAT_GNN_DATA_DP, "cat_topk.txt")
    SOL_GNN_TOPK_FP = os.path.join(_SOL_GNN_DATA_DP, "sol_topk.txt")

    CAT_GNN_TOPK_ACC_FP = os.path.join(_CAT_GNN_DATA_DP, "cat_topk_acc.tsv")
    SOL_GNN_TOPK_ACC_FP = os.path.join(_SOL_GNN_DATA_DP, "sol_topk_acc.tsv")

    CAT_GNN_MODEL_DP = os.path.join(_CAT_GNN_DATA_DP, "models")
    SOL_GNN_MODEL_DP = os.path.join(_SOL_GNN_DATA_DP, "models")

    TENSOR_BOARD_DP = os.path.join(_SS_APP_DP, "logs")

    CLEAN_RXN_ERROR_FP = os.path.join(_SS_APP_DP, "clean_rxn_error.tsv")

    dirs = [_SS_APP_DP, _CAT_GNN_DATA_DP, _SOL_GNN_DATA_DP,
            CAT_GNN_SINGLE_DATA_DP, SOL_GNN_SINGLE_DATA_DP,
            TRAIN_CAT_GNN_BATCH_DATA_DP, TRAIN_SOL_GNN_BATCH_DATA_DP,
            TEST_CAT_GNN_BATCH_DATA_DP, TEST_SOL_GNN_BATCH_DATA_DP,
            CAT_GNN_MODEL_DP, SOL_GNN_MODEL_DP,
            TENSOR_BOARD_DP]


class BuildingBlocksPath(Path):

    _BB_DP = os.path.join(_PROJ_PATH, "BuildingBlocks")

    BBM_FP = os.path.join(_BB_DP, "rxndb_public_building_block_mols.tsv")
    BBM_NON_STEREO_FP = os.path.join(_BB_DP, "non_stereo_rxndb_public_building_block_mols.tsv")

    SM_FP = os.path.join(_BB_DP, "rxndb_public_for_sale_mols.tsv")
    SM_NON_STEREO_FP = os.path.join(_BB_DP, "non_stereo_rxndb_public_for_sale_mols.tsv")

    RM_FP = os.path.join(_BB_DP, "rxndb_public_reactant_mols.tsv")
    RM_NON_STEREO_FP = os.path.join(_BB_DP, "non_stereo_rxndb_public_reactant_mols.tsv")

    dirs = []


class RxnDeterminerPath(Path):

    RXN_WITH_SPECTATORS_FP = SpectatorsSelectorPath.RXN_WITH_SPECTATORS_FP
    CAT_GNN_DATA_INFO_FP = SpectatorsSelectorPath.CAT_GNN_DATA_INFO_FP
    SOL_GNN_DATA_INFO_FP = SpectatorsSelectorPath.SOL_GNN_DATA_INFO_FP

    _RD_DP = os.path.join(_PROJ_PATH, "RxnDeterminer")

    RIGHT_RXN_FP = os.path.join(_RD_DP, "right_rxn.tsv")
    RIGHT_AND_YIELD_RXN_FP = os.path.join(_RD_DP, "right_and_yield_rxn.tsv")
    RIGHT_AND_HIGH_YIELD_RXN_FP = os.path.join(_RD_DP, "right_and_high_yield_rxn.tsv")

    HIGH_YIELD_SPEC_RXN_FP = os.path.join(_RD_DP, "high_yield_spec_rxn.tsv")

    WRONG_RXN_FP = os.path.join(_RD_DP, "wrong_rxn.tsv")
    WRONG_RXN_FPS = os.path.join(_RD_DP, "wrong_rxn_{0}.tsv")

    RXN_RIGHT_SINGLE_DATA_DP = os.path.join(_RD_DP, "rxn_right_single_data")
    RXN_WRONG_SINGLE_DATA_DP = os.path.join(_RD_DP, "rxn_wrong_single_data")
    RXN_TRAIN_SINGLE_FNS_FP = os.path.join(_RD_DP, "rxn_train_single_fns.txt")
    RXN_TEST_SINGLE_FNS_FP = os.path.join(_RD_DP, "rxn_test_single_fns.txt")
    RXN_TRAIN_BATCH_GNN_DATA_DP = os.path.join(_RD_DP, "rxn_train_batch_gnn_data")
    RXN_TEST_BATCH_GNN_DATA_DP = os.path.join(_RD_DP, "rxn_test_batch_gnn_data")

    MODEL_CAT_SOL_DP = os.path.join(_RD_DP, 'model_cat_sol')
    MODEL_NO_CAT_SOL_DP = os.path.join(_RD_DP, 'model_no_cat_sol')

    dirs = [_RD_DP, RXN_TRAIN_BATCH_GNN_DATA_DP, RXN_TEST_BATCH_GNN_DATA_DP,
            RXN_RIGHT_SINGLE_DATA_DP, RXN_WRONG_SINGLE_DATA_DP, MODEL_CAT_SOL_DP,
            MODEL_NO_CAT_SOL_DP]


BaseTempPath.init()
CompTempPath.init()
SpectatorsSelectorPath.init()
RxnDeterminerPath.init()


if __name__ == '__main__':
    pass
