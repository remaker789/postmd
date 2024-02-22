#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''This script is to process the files generated from log file(log.lammps).
@Author  :   Shusong Zhang
@Email   :   sszhang@mail.nwpu.edu.cn, zhangshusong789@gmail.com
@File    :   FileOperation.py
@Time    :   2023/05/05 10:46:11
'''

import re, os
from pathlib import Path
from ..utils import judge_file


class LogFile:
    def __init__(self) -> None:
        self.lammps_version=None
        self.path=None
        
     
        
    def set_path(self, path=None):
        """set the path of log file(log.lammps) and judge the validity of the path

        Args:
            path (str, optional): the path to the log file. Defaults to None.
        """        
        self.path = path
        judge_file(self.path)
        print(f"You are processing the file: '{self.path}'")

 
        
        
    def get_lammps_version(self):
        """get the version of lammps used.
        """     
        judge_file(self.path)   
        with open(self.path, encoding='utf-8') as f:
            self.lammps_version = f.readline()
            print(f"The version of LAMMPS: {self.lammps_version}")
    

    def extract_thermodata(self,*, block_start = "Per MPI rank", block_end = "Loop time of", path=None, output="extracted-log.lammps"):        
        """extract the thermo data from log file

        Args:
            block_start (str, optional): the start string of the block of thermo data. Defaults to "Per MPI rank".
            block_end (str, optional): the end string of the block of thermo data. Defaults to "Loop time of".
            paths (str, optional): the path to the log file. Default to None, which means use self.path
            
        .. note::
            1. block_start='Step' and block_end="Loop time of" usually works for most case in LAMMPS (29 Oct 2020). 
            If the block changed with lammps, you can change it with keyword arguments.

            2. If the block of thermo data contain the WARNING, for example, sub-domain..., this function can not omit the WARNING.
            Maybe a good choice is to add a judgement depend on the percent of characters exceeding 50% with Regular Expression!

        """
        
        path = path if path else self.path
        self.set_path(path)
        output_path = os.path.join(os.path.dirname(path), output)
        with open(self.path, encoding='utf-8') as logfile:
            lines_block_start=[] # the list storing the line number of block start
            lines_block_end=[]   # the list storing the line number of block end
            num_block=0 # count the number of blocks
            line_block_start=0
            line_block_end=0
            write=False
            text=""
            for i, line in enumerate(logfile.readlines()):
                # todo 看看能不能再精简代码结构。。。现在还是有点麻烦
                if line.startswith(block_end):
                    write=False
                    line_block_end=i-1
                    lines_block_end.append(line_block_end+1) # the line in python starts from 0

                if not line.startswith(block_start):
                    if  i == line_block_start and len(lines_block_end) != 0: # avoid output the header second time.
                        continue
                    
                    if write:
                        text=text+line
                    else:
                        continue
                else:
                    write=True
                    line_block_start=i+1
                    lines_block_start.append(line_block_start+1) # the line in python starts from 0

                
            output_path= os.path.join(os.path.dirname(self.path), output)
            with open(output_path, 'w+') as output:
                output.write(text)
                
        print(f"The line number of block start: {lines_block_start}")   
        print(f"The line number of block end:   {lines_block_end}")   
        print("---------- Extracting logfile succeed! -----------")
        print(f"The extracted file: '{output_path}'")
