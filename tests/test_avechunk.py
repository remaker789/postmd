#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@Author  :   Shusong Zhang
@Email   :   sszhang@mail.nwpu.edu.cn, zhangshusong789@gmail.com
@File    :   test_avechunk.py
@Time    :   2024/03/03
'''

import os
os.chdir(os.path.dirname(__file__))

from postmd.avechunk import Bin1d


filepath = r"data/ave_chunk_bin1d.dat"

bin1d = Bin1d(filepath,dim="z")
bin1d.read_file()

# We stored the blocks in a list
print(f"The type of bin1d.blocks is {type(bin1d._blocks)}")

# you can use len(bin1d) to get the number of blocks
print(f"{bin1d.path} have {len(bin1d)} blocks")

# you can use bin1d[index] to access index-th block
# in each block, we store the information with dict.
print("The first block is:")
print(bin1d[0])

# you can use key to access the data in each block
print("The data in the first block is:")
print(bin1d[0]["data"])




