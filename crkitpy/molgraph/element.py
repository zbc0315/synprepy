# -*- coding: utf-8 -*-
# @Time     : 2020/4/20 22:22
# @Author   : zbc
# @Email    : zhangbc0315@outlook.com
# @FileName : element.py
# @Software : PyCharm


class Element:

    elements = ['h', 'he',
                'li', 'be', 'b', 'c', 'n', 'o', 'f', 'ne',
                'na', 'mg', 'al', 'si', 'p', 's', 'cl', 'ar',
                'k', 'ca', 'sc', 'ti', 'v', 'cr', 'mn', 'fe', 'co', 'ni', 'cu', 'zn', 'ga', 'ge', 'as', 'se', 'br', 'kr',
                'rb', 'sr', 'y', 'zr', 'nb', 'mo', 'tc', 'ru', 'rh', 'pd', 'ag', 'cd', 'in', 'sn', 'sb', 'te', 'i', 'xe',
                'cs', 'ba', 'la', 'ce', 'pr', 'nd', 'pm', 'sm', 'eu', 'gd', 'tb', 'dy', 'ho', 'er', 'tm', 'yb', 'lu', 'hf', 'ta', 'w', 're', 'os', 'ir', 'pt', 'au', 'hg', 'tl', 'pb', 'bi', 'po', 'at', 'rn',
                'fr', 'ra', 'ac', 'th', 'pa', 'u', 'np', 'pu', 'am', 'cm', 'bk', 'cf', 'es', 'fm', 'md', 'no', 'lr', 'rf', 'db', 'sg', 'bh', 'hs', 'mt', 'ds', 'rg', 'cn', 'nh', 'fl', 'mc', 'lv', 'ts', 'og', 'uun', 'uuu', 'uub', 'uut', 'uuq', 'uup', 'uuh', 'uus', 'uuo']

    @classmethod
    def get_atomic_num_by_symbol(cls, symbol):
        return cls.elements.index(symbol.lower()) + 1
