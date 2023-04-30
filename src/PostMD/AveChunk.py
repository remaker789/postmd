import re

import numpy as np
import pandas as pd

'''
This script is to deal with the data using comput chunk/atom and fix ave/chunk.
now I just right the Bin1d Class
- This script can only be used for constant chunk numbers. using a iter to save memory.
- For changed chunk numbers, I may used skip_rows..
'''

class Bin1d:
    def __init__(self, dim:str, delta:float):
        """initialize the bin1s

        Args:
            dim (str): the dimension of chunk. 'x','y' or 'z'
            delta (float): the chunk size in Angstrom
        """        
        self.dim = dim
        self.delta = delta
        
        
    def get_names(self, num_line:int)->list: 
        """get the names used in dataframe when read files

        Args:
            num_line (int): the line including the column names of data

        Returns:
            list: a list of the column names of data
        """        
        with open(self.filepath) as f:
            for i, line in enumerate(f.readlines(),1):
                if i == num_line:
                    return re.split('\s+', line.strip())[1:]
                
                
    def read_file(self,
                  filepath:str, 
                  num_line_time:int=2, 
                  names:list=None, 
                  num_line_column:int=3
                  , **kwargs)->pd.DataFrame:
        """read file and return DataFrame. 
            There are some default setting in read files:
            - comment="#" due to the nature of LAMMPS
            - sep="\s+" to match all space between data
            - **kwargs received parameters of read_csv and read_excel, depending on cases. 

        Args:
            filepath (str): path to the file
            names (list, optional): list of columns names of data. 
                                    Defaults to None, meaning using the content in num_line of file.
            num_line_time (int, optional): "# Timestep Number-of-chunks Total-count" For bin/1d. Default=2
            num_line_column (int, optional): when names=None, the content in line=num_line
                                      will be used as the column names. Defaults to 3.

        Returns:
            pd.DataFrame: DataFrame object read from the file.
        """
        time_line = self.get_names(num_line_time)       
        names = names if names else self.get_names(num_line_column) 

        df = pd.read_csv(filepath, names =names, sep = '\s+', comment = "#", **kwargs)
        return df   
    
    
    # chunk data的文件总体结构为：不同时间段内进行的统计，分为一个个block，因此我们将block的数目设置为len.
    def __len__(self):
    # 将总共进行了多少次chunk统计设置为len()    
        pass
    def __getitem__(self):
        # 索引第n次的chunk结果
        pass