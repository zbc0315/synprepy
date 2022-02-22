# -*- coding: utf-8 -*-
# @Time     : 2021/4/16 14:36
# @Author   : zbc
# @Email    : zhangbc0315@outlook.com
# @FileName : note.py
# @Software : PyCharm

import os


class Note:

    @classmethod
    def note(cls, note_fp: str, line: str) -> None:
        mode = 'a'
        if not os.path.exists(note_fp):
            mode = 'w'
        with open(note_fp, mode, encoding='utf-8')as f:
            f.write(line)
            f.write('\n')


if __name__ == '__main__':
    pass
