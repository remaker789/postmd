#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@Author  :   Shusong Zhang
@Email   :   sszhang@mail.nwpu.edu.cn, zhangshusong789@gmail.com
@File    :   test_greenkubo.py
@Time    :   2023/04/27 15:34:34
'''

import matplotlib.pyplot as plt

raw_path=r"./data/friction.txt"

from PostMD.GreenKubo import *

test1=GreenKubo()
test1.read_file(raw_path)
test1.cal_acf(data_type="raw",col=3,nlag=2000)
plt.plot(test1.nlag, test1.acf,label="from raw data")

ave_path=r"./data/friction_acf_overwriting.txt"
test2=GreenKubo()
test2.read_file(ave_path,name_line=3,skiprows=4)
test2.cal_acf(data_type='acf',col=5)
plt.plot(test2.nlag,test2.acf, label="from fix ave/correlate")

plt.xlabel("nlag")
plt.ylabel("acf")
plt.legend()


## ---------- integrate acf ------------
plt.figure()

test1.integrate_acf()
plt.plot(test1.nlag, test1.int_acf, label="simpson method from raw data")

test1.integrate_acf(method="trap")
plt.plot(test1.nlag, test1.int_acf, label="trap method from raw data")

test2.integrate_acf()
plt.plot(test2.nlag, test2.int_acf, label="simpson method from fix ave/correalte")

test2.integrate_acf(method="trap")
plt.plot(test2.nlag, test2.int_acf, label="trap method from fix ave/correalte")

plt.xlabel("nlag")
plt.ylabel("integration of acf")
plt.legend()
plt.show()

## the acf in my work
path_mywork=r"./data/friction_acf_overwriting_mywork.txt"
mywork=GreenKubo()
mywork.read_file(path=path_mywork, name_line=3, skiprows=4)
mywork.cal_acf(data_type="acf",col=5,nlag_col=1)
plt.figure()
plt.plot(mywork.nlag,mywork.acf)

plt.xlabel("nlag")
plt.ylabel("acf")
plt.legend()

plt.figure()
mywork.integrate_acf(method="simpson")
plt.plot(mywork.nlag, mywork.int_acf, label="simpson method")
mywork.integrate_acf(method="trap")
plt.plot(mywork.nlag, mywork.int_acf, label="trap method")
plt.xlabel("nlag")
plt.ylabel("integration of acf")
plt.legend()
plt.show()
