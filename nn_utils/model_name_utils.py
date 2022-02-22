#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Time    : 2021/8/24 14:33
# @Author  : zhangbc0315@outlook.com
# @File    : model_name_utils.py
# @Software: PyCharm


class ModelNameUtils:

    """ 机器学习模型的命名与名称解析

    """

    @classmethod
    def get_model_name(cls, idx: int, loss: float, acc: float):
        return f"{idx}_{loss}_{acc}_.pt"

    @classmethod
    def parse_model_name(cls, model_name: str):
        split_name = model_name.split('_')
        res = {'idx': int(split_name[0]),
               'loss': float(split_name[1]),
               'acc': float(split_name[2])}
        return res

    @classmethod
    def get_fn_from_fp(cls, fp: str):
        return fp.split('/')[-1].split('\\')[-1]

    @classmethod
    def parse_model_file_path(cls, model_fp: str):
        return cls.parse_model_name(cls.get_fn_from_fp(model_fp))


if __name__ == "__main__":
    print(ModelNameUtils.get_fn_from_fp("C:\\asd/asd/asf\\asd.fp"))
