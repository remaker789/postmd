#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@Author  :   Shusong Zhang
@Email   :   sszhang@mail.nwpu.edu.cn, zhangshusong789@gmail.com
@File    :   test_greenkubo.py
@Time    :   2023/04/27 15:34:34
'''
import numpy as np
import matplotlib.pyplot as plt
import os
import scipy.constants as C
os.chdir(os.path.dirname(__file__))

from postmd.avetime import GreenKubo

e=C.e

# the acf data, included two forms: 
# one is the raw data of friction force at each step of CNT, in file "friction.txt"
# another is the autocorrelated data of friction force by 'fix ave/correlate', in file "friction_acf_overwriting.txt"
raw_data_path=r"./data/friction.txt"
acf_data_path=r"./data/friction_acf_overwriting.txt"




## ================ For raw data ===================
unit_trans = e/1e-10  # friction force的单位是eV/Ang
mywork=GreenKubo(raw_data_path)
mywork.read_file()
mywork.calc_acf(data_type="raw",col=3,nlag=2000,unit_trans=unit_trans)
plt.figure()
plt.plot(mywork.nlag,mywork.acf)
plt.title("computed from raw data")
plt.xlabel("nlag")
plt.ylabel("acf")
plt.legend()
# plt.savefig("./output/green_kubo_raw_data_acf.png")


A = 2*np.pi*2.335e-10*200e-10
kB = C.Boltzmann
T = 298
plt.figure()
mywork.integrate_acf(method="simpson")
plt.plot(mywork.nlag, mywork.int_acf/(A*kB*T), label="simpson method")
mywork.integrate_acf(method="trap")
plt.plot(mywork.nlag, mywork.int_acf/(A*kB*T), label="trap method")
# print( (mywork.int_acf/(A*kB*T)).max())
plt.xlabel("nlag")
plt.ylabel("integration of acf")
plt.title("computed from raw data")
plt.legend()
# plt.savefig("./output/green_kubo_raw_data_int_acf.png")



## ================ For autocorrelation data ==================
unit_trans = e/1e-10  # friction force的单位是eV/Ang
mywork=GreenKubo(path=acf_data_path)
mywork.read_file(header_line=3, skiprows=4)
mywork.calc_acf(data_type="acf",col=5,nlag_col=1,unit_trans=unit_trans)
plt.figure()
plt.plot(mywork.nlag,mywork.acf)
plt.title("computed from acf data")
plt.xlabel("nlag")
plt.ylabel("acf")
plt.legend()
# plt.savefig("./output/green_kubo_acf_data_acf.png")

A = 2*np.pi*2.335e-10*200e-10
kB = C.Boltzmann
T = 298
plt.figure()
mywork.integrate_acf(method="simpson")
plt.plot(mywork.nlag, mywork.int_acf/(A*kB*T), label="simpson method")
mywork.integrate_acf(method="trap")
plt.plot(mywork.nlag, mywork.int_acf/(A*kB*T), label="trap method")
# print( (mywork.int_acf/(A*kB*T)).max())
plt.xlabel("nlag")
plt.ylabel("integration of acf")
plt.title("computed from acf data")
plt.legend()
# plt.savefig("./output/green_kubo_acf_data_int_acf.png")



plt.show()


# You will find the output is same for both acf and raw data, because these two files are generated from the same MD simulation.
# This means you can use "fix ave/correlate" to process acf on the fly, or calculate acf by post-processing the raw data.
