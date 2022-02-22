# 反应判断器

### 目的
根据反应物和产物，判断反应是否合理

### 一些概念
SMILES：通过一串字符，表示一个分子的结构，例如：“C1CCCCC1”表示环己烷

INCHI：和SMILES作用类似，只是字符的规则不同，例如：“InChI=1S/C6H12/c1-2-4-6-5-3-1/h1-6H2”表示环己烷

不需要我们理解这些字符串的构成规则，有一些软件可以解析或导出，例如rdkit

```python
from rdkit.Chem import AllChem

mol = AllChem.MolFromSmiles("C1CCCCC1")
```

### 数据

#### 合理的反应数据

rxndb_public_clean_duplicate_rxn.tsv

表中的rxn_smi即为反应的SMILES，反应物的SMILES是由化合物的SMILES组合而成，规则是：

假如反应是A+B->C,

则反应的SMILES为“A.B>>C”

假如反应是A+B->C(溶剂/催化剂为D和F)

则反应的SMILES为“A.B>D.F>C”

这些规则也不需要掌握，同样可以使用rdkit解析或导出

```python
from rdkit.Chem import AllChem

rxn = AllChem.ReactionFromSmarts("C.CC>>CCC")
```

#### 不合理的反应数据

不合理的反应数据存储在两个文件中：

rxndb_wrong_rxn_wrong_rxn.tsv

rxn_code给出反应，例如14973.157986>>49684，这些数字是mid，但未给出反应物产物的INCHI或SMILES；

rxndb_wrong_rxn_wrong_rxn_mols.tsv

给出了化合物的INCHI

### 任务

使用全连接神经网络，预测反应是否合理

提示：可使用rdkit将化合物及化学反应转化为分子指纹，作为模型输入

