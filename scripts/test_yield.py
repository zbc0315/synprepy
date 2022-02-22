#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Time    : 2021/9/16 10:06
# @Author  : zhangbc0315@outlook.com
# @File    : test_yield.py
# @Software: PyCharm


def mcts():
    r = None
    print('开始')
    while True:
        n = yield f'收到了n{r}'
        r = n
        if n == 2:
            print('收到停止信号')
            yield f'停止后的消息'
            break


def test():
    m = mcts()
    m.send(None)
    for n in range(3):
        res = m.send(n)
        print(res)
    m.close()


if __name__ == "__main__":
    test()
