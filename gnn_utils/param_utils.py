#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Time    : 2021/8/18 13:24
# @Author  : zhangbc0315@outlook.com
# @File    : param_utils.py
# @Software: PyCharm


class ParamUtils:

    @classmethod
    def adjust_lr(cls, lr_start, lr_end, optimizer, epoch_num, lr_change_epoch_num=30):
        lr = lr_start * (0.1 ** (epoch_num // lr_change_epoch_num))
        lr = max(lr, lr_end)
        print(f"epoch: {epoch_num}, change lr to : {lr}")
        for param_group in optimizer.param_groups:
            param_group['lr'] = lr


if __name__ == "__main__":
    pass
