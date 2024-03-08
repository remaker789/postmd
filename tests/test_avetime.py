#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@Author  :   Shusong Zhang
@Email   :   sszhang@mail.nwpu.edu.cn, zhangshusong789@gmail.com
@File    :   tset.py
@Time    :   2023/04/27 15:16:21
'''
import os 
os.chdir(os.path.dirname(__file__))

from postmd.avetime import AveTime

test = AveTime(r"./data/ave_time.dat")
df = test.read_file()
print(df.head())