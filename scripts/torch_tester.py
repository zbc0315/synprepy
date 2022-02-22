# -*- coding: utf-8 -*-
# @Time     : 2021/5/18 10:58
# @Author   : zbc
# @Email    : zhangbc0315@outlook.com
# @FileName : torch_tester.py
# @Software : PyCharm

import torch


class TorchTester:

    @classmethod
    def test_sort(cls):
        d = torch.tensor([[7, 2, 3], [0.1, 0.5, 0.2]], dtype=torch.float)
        sort_d = torch.sort(d, dim=1, descending=True)
        sort_d_idx = sort_d[1]

        y = torch.tensor([[0, 0.9, 0], [0, 0, 0.9]], dtype=torch.float)
        max_y = y.max(dim=1)
        max_y_idx = y.max(dim=1)[1]
        print("end")

    @classmethod
    def test_len(cls):
        s = torch.tensor([0, 0, 1], dtype=torch.float)
        print(s.size().numel())


if __name__ == '__main__':
    TorchTester.test_sort()
    # TorchTester.test_len()
