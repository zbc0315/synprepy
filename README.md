
# 化学反应逆向预测

## Python 版本

python3.6

## 安装依赖
* RDKit
```bash
# 安装rdkit
conda install -c conda-forge rdkit
# 安装pytorch
请按照官网安装: https://pytorch.org/
# 安装pyg
请按照官网安装: https://pytorch-geometric.readthedocs.io/en/latest/notes/installation.html
# 安装其他依赖
cd synprepy # 进入项目文件夹
pip install -r requirements.txt

```

## 准备数据

1. 准备具有如下反应数据的tsv文件

| rid | rxn_smi |
| --- | ------- |
| 0 | [OH:2][N:3]1[N:4]=[N:5][C:6]2=[C:7]1[CH:8]=[CH:9][CH:10]=[CH:11]2>>[OH:2][N:3]1[N:4]=[N:5][C:6]2=[C:7]1[CH:8]=[CH:9][CH:10]=[CH:11]2 |
| 1 | [O:13]1[CH2:14][CH2:15][CH2:16][CH2:17]1>>[O:13]([CH2:14][C@@H:15]1[O:13][C@H:17]2[C@@H:15]([CH2:14][CH2:16]1)[CH2:16]2)[CH3:17] |


## 填写配置文件（config.json）

```json
{
    "root": "存储反应数据和反应模板数据的文件夹路径",
    "rxn_data_tsv_file_name": "反应数据文件名（该文件必须处于上个设置的root文件夹内）",
    "rid_with_rxn_template_tsv_file_name": "保持默认即可",
    "rxn_centralized_template_tsv_file_name": "保持默认即可",
    "rxn_extended_template_tsv_file_name": "保持默认即可",
    "rid_with_tid_tsv_file_name": "保持默认即可",
    "rxn_centralized_template_selector_config": {
        "root": "存储中心反应模板选择器的相关文件的文件夹",
        "train_rids_file_name": "保持默认即可",
        "test_rids_file_name": "保持默认即可",
        "train_temp_dir_name": "保持默认即可",
        "test_temp_dir_name": "保持默认即可",
        "min_num_covered_rxns_by_rxn_centralized_template": "过滤模板时要求模板的最少覆盖反应数",
        "filter_tids_file_name": "保持默认即可",
        "device": "cpu或者cuda",
        "batch_size": 1000,
        "epoch_num": 200,
        "lr_start": 0.0005,
        "lr_end": 1e-7
    },

  "min_num_covered_rxns_by_rxn_extended_template": 50
}
```

## 单步反应预测

1. 抽取反应模板：

运行 rxn_template_prepare.py

2. 训练单步反应模板选择器：

运行 rxn_template_selector/rct_train_test.py
