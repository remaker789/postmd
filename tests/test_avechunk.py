#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@Author  :   Shusong Zhang
@Email   :   sszhang@mail.nwpu.edu.cn, zhangshusong789@gmail.com
@File    :   test_avechunk.py
@Time    :   2023/08/24 11:09:42
'''

import os
from PostMD.AveChunk import *
os.chdir(os.path.dirname(__file__))

filepath = r"data/ave_chunk_bin1d.dat"

bin1d = Bin1d(dim="z")
bin1d.read_file(filepath=filepath)
print(type(bin1d.content))
print(bin1d.content[0])
print(len(bin1d))
print(bin1d[1])
print(bin1d[2]["data"])

