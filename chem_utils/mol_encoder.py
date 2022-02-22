# -*- coding: utf-8 -*-
# @Time     : 2021/4/12 10:41
# @Author   : zbc
# @Email    : zhangbc0315@outlook.com
# @FileName : mol_encoder.py
# @Software : PyCharm

from rdkit.Chem import AllChem


class MolEncoder:

    @classmethod
    def _ten_to_nine(cls, num: int) -> [int]:
        """ 十进制转化为九进制

        :param num:
        :return:
        """
        res = []
        while True:
            n = num % 9
            res.append(int(n))
            rest = (num - n) / 9
            if rest == 0:
                break
            num = rest
        res.reverse()
        return res

    @classmethod
    def _remove_0(cls, num: int) -> str:
        """ 去除数字中的0

        :param num:
        :return:
        """
        nine_nums = cls._ten_to_nine(num)
        nine_nums_no_0_str = [str(n + 1) for n in nine_nums]
        return ''.join(nine_nums_no_0_str)

    @classmethod
    def mol_to_code_by_smiles(cls, smiles: str):
        mol = AllChem.MolFromSmiles(smiles)
        atom_codes = []
        for atom in mol.GetAtoms():
            atom_code_nums = [cls._remove_0(atom.GetAtomicNum()),
                              cls._remove_0(atom.GetTotalNumHs()),
                              cls._remove_0(atom.GetFormalCharge()),
                              cls._remove_0(int(atom.GetChiralTag()))]
            atom_code = '0'.join([str(c) for c in atom_code_nums])
            atom_codes.append(atom_code)
        bond_codes = []
        for bond in mol.GetBonds():
            bond_code_nums = [cls._remove_0(int(bond.GetBondType())),
                              cls._remove_0(int(bond.GetBondDir())),
                              cls._remove_0(bond.GetBeginAtomIdx()),
                              cls._remove_0(bond.GetEndAtomIdx())]
            bond_code = '0'.join([str(c) for c in bond_code_nums])
            bond_codes.append(bond_code)
        return '00'.join(atom_codes) + '000' + '00'.join(bond_codes)


if __name__ == '__main__':
    r = MolEncoder.mol_to_code_by_smiles('c1ccccc1')
    print(r)
