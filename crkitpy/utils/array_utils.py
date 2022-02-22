# -*- coding: utf-8 -*-
# @Time     : 2020/4/22 10:23
# @Author   : zbc
# @Email    : zhangbc0315@outlook.com
# @FileName : array_utils.py
# @Software : PyCharm

from typing import Dict, List


class ArrayUtils:

    @classmethod
    def select_random_array_in_matrix(cls, matrix_dict: Dict[object, List[object]]) -> Dict[object, object]:
        """ 处理所有value都为list的dict， 从每一个list形式的value中抽取一个值作为当前key的新value，最后形成的新dict有可能会出现相同的value

        :param matrix_dict:
        :return:
        """
        array_dict = {}
        selected_values = []
        for key, values in matrix_dict.items():
            has_choice = False
            for value in values:
                if value not in selected_values:
                    array_dict[key] = value
                    selected_values.append(value)
                    has_choice = True
            if not has_choice:
                array_dict[key] = values[-1]
        return array_dict

    @classmethod
    def select_unique_array_in_matrix(cls, matrix_dict: Dict[object, List[object]]) -> Dict[object, object]:
        """ 处理所有value都为list的dict， 从每一个list形式的value中抽取一个值作为当前key的新value，确保最后形成的新dict中每一个key对应的value都是不同的

        :param matrix_dict: 所有value都为list的dict, 例如: matrix_dict = {'k1': ['a', 'b'],
                                                                        'k2': ['b', 'c'],
                                                                        'k3': ['c']}
        :return: 处理之后的新dict，例如: array_dict = {‘k1’: 'a',
                                                    'k2': 'b',
                                                    'k3': 'c'}
        """
        all_values = []
        for key, values in matrix_dict.items():
            for value in values:
                if value not in all_values:
                    all_values.append(value)
        if len(all_values) < len(matrix_dict.keys()):
            return {}

        matrix_dict = cls._sort_dict_by_value_length(matrix_dict)
        selected_values = []
        cls._select_unique_array(matrix_dict, 0, selected_values)
        array_dict = {}
        for n, key in enumerate(matrix_dict.keys()):
            array_dict[key] = selected_values[n]
        return array_dict

    @classmethod
    def _select_unique_array(cls, matrix_dict: Dict[object, List[object]], current_idx: int, selected_values: List[object]):
        key = list(matrix_dict.keys())[current_idx]
        values = matrix_dict[key]
        for selected_value in selected_values:
            if selected_value in values:
                values.remove(selected_value)
        if len(values) == 0:
            return False
        for value in values:
            selected_values.append(value)
            if current_idx == len(matrix_dict.keys()) - 1:
                return True
            if cls._select_unique_array(matrix_dict, current_idx + 1, selected_values):
                return True
            else:
                selected_values = selected_values[:current_idx]
        return False

    @staticmethod
    def _sort_dict_by_value_length(dic: Dict[object, List[object]]) -> Dict[object, List[object]]:
        """ 依照value的长度对dict进行排序

        :param dic: 排序前的字典，例如 dic = {'a': ['a1', 'a2'],
                                           'b': ['b1', 'b2', 'b3'],
                                           'c': ['c1']}
        :return: 排序后的字典，例如 new_dic = {'c': ['c1'],
                                            'a': ['a1', 'a2'],
                                            'b': ['b1', 'b2']}
        """
        items = sorted(dic.items(), key=lambda x: len(x[1]))
        new_dic = {}
        for item in items:
            new_dic[item[0]] = item[1]
        return new_dic
