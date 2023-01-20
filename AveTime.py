import numpy as np
import pandas as pd
from pathlib import Path
import re

class AveTime:
    def __init__(self, filepath=None, timestep = 1):  # timestep = 1fs
        self.timestep = timestep
        if filepath is None:
            raise ValueError("filepath can not be empty")
        if not Path(filepath).is_file():
            raise ValueError(f"{filepath} is not a file")
        self.filepath = filepath

        
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
            
            
            
            
    def read_file(self,filepath:str, names:list=None, num_line:int = 2, **kwargs)->pd.DataFrame:
        """read file and return DataFrame. 
            There are some default setting in read files:
            - comment="#" due to the nature of LAMMPS
            - sep="\s+" to match all space between data
            - **kwargs received parameters of read_csv and read_excel, depending on cases. 

        Args:
            filepath (str): path to the file
            names (list, optional): list of columns names of data. 
                                    Defaults to None, meaning using the content in num_line of file.
            num_line (int, optional): when names=None, the content in line=num_line
                                      will be used as the column names. Defaults to 2.

        Returns:
            pd.DataFrame: DataFrame object read from the file.
        """        
        names = names if names else self.get_names(num_line)
        # if filepath.endswith(('xls', 'xlsx')): #产生的数据基本上不可能是excel
        #     df = pd.read_excel(filepath, names =names, comment = "#", **kwargs) 
        # else:
        #     df = pd.read_csv(filepath, names =names, sep = '\s+', comment = "#", **kwargs)
        
        df = pd.read_csv(filepath, names =names, sep = '\s+', comment = "#", **kwargs)
        return df



path = r"D:\OneDrive\Github\MyRepos\PostMD\data\flux.dat"

test = AveTime(path)
print(test.get_names(2))
df = test.read_file(path)
print(df.head(5))