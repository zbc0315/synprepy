# 反应/分子预测

## 安装依赖
```bash
conda install -c rdkit rdkit
conda install -c openbabel openbabel
pip install -r requirements.txt
```
安装pytorch，参考官网，选择版本(1.9)
https://pytorch.org/get-started/locally/

安装pytorch_geometric，参考文档，选择版本(注意要和pytorch版本匹配)
https://pytorch-geometric.readthedocs.io/en/latest/notes/installation.html


预测催化剂，得cid，长度为催化剂数目：325

催化剂，不使用催化剂:idx=0，使用未知催化剂:idx=num_cats，使用催化剂时，cid+1